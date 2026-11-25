"""Measure population in |2> after one X01 pulse on a Qolab Start transmon."""

import argparse
import json
from pathlib import Path

import numpy as np


STATES = ("g", "e", "f", "test")
ANALYSES = ("map2", "none", "map1")
REF_OP = "p2_reference_01"
EF_OP = "p2_reference_12"


def build_program(machine, qubit, args):
    from qm.qua import FUNCTIONS, align, declare, declare_output_stream, for_, program, save, stream_processing
    from quam.components.pulses import DragCosinePulse

    if args.test_op not in qubit.xy.operations or args.test_op in (REF_OP, EF_OP):
        raise ValueError("Select an existing test operation on the qubit XY channel.")

    # Attach two independently calibrated cosine pulses to the in-memory state.
    # Alpha=0 gives plain cosine envelopes for these deliberately selective pulses.
    for name, amplitude, duration in (
        (REF_OP, args.ref_amp, args.ref_ns),
        (EF_OP, args.ef_amp, args.ef_ns),
    ):
        if name in qubit.xy.operations:
            raise ValueError(f"Reserved operation already exists: {name}")
        qubit.xy.operations[name] = DragCosinePulse(
            length=duration, amplitude=amplitude, axis_angle=0.0,
            alpha=0.0, anharmonicity=0.0, digital_marker="ON",
        )

    config = machine.generate_config()
    ge_if = round(qubit.xy.intermediate_frequency)
    ef_if = round(ge_if + args.f12_hz - qubit.xy.RF_frequency)
    if ge_if == ef_if:
        raise ValueError("The GE and EF drive frequencies must differ.")
    analyses = ANALYSES if args.full else ANALYSES[:1]

    def play_ef():
        qubit.xy.update_frequency(ef_if)
        qubit.xy.play(EF_OP)
        qubit.xy.update_frequency(ge_if)

    with program() as experiment:
        shot = declare(int)
        state = declare(int)
        streams = {(a, s): declare_output_stream() for a in analyses for s in STATES}
        machine.initialize_qpu(target=qubit)
        align()
        with for_(shot, 0, shot < args.blocks * args.shots_per_block, shot + 1):
            for analysis in analyses:
                for preparation in STATES:
                    qubit.align()
                    qubit.reset("thermal")
                    qubit.align()
                    if preparation in ("e", "f"):
                        qubit.xy.play(REF_OP)
                    if preparation == "f":
                        play_ef()
                    if preparation == "test":
                        qubit.xy.play(args.test_op)
                    if analysis == "map2":
                        play_ef()
                    if analysis in ("map1", "map2"):
                        qubit.xy.play(REF_OP)
                    qubit.align()
                    qubit.readout_state(state)
                    save(state, streams[analysis, preparation])
        with stream_processing():
            for (analysis, preparation), stream in streams.items():
                stream.buffer(args.shots_per_block).map(FUNCTIONS.average()).save_all(f"{analysis}_{preparation}")
    return experiment, config, analyses, ge_if, ef_if


def estimate(bright, background_tolerance, draws=4000):
    """Return a signed estimate and a block-bootstrap interval in probability units."""
    dark = 1.0 - np.asarray(bright, dtype=float)
    if (dark.ndim != 3 or dark.shape[0] < 8 or dark.shape[1:] not in ((1, 4), (3, 4))
            or not np.isfinite(dark).all() or np.any((dark < 0) | (dark > 1))):
        raise ValueError("Expected at least eight blocks of valid bright fractions.")
    rng = np.random.default_rng(763)
    indices = rng.integers(len(dark), size=(draws, len(dark)))
    samples = dark[indices].mean(axis=1)
    mean = dark.mean(axis=0)
    result = {"dark_fractions": mean.tolist()}

    if dark.shape[1] == 1:
        g, e, f, test = mean[0]
        bg = (g + e) / 2
        g_b, e_b, f_b, t_b = samples[:, 0].T
        delta = float(np.quantile(abs(g_b - e_b), 0.95))
        contrast = float(np.quantile(f_b - np.maximum(g_b, e_b), 0.025))
        result.update(background_difference_95=delta, mapping_contrast_lower_95=contrast)
        if contrast <= 0:
            return result | {"status": "mapping_contrast_unresolved"}
        if delta > background_tolerance:
            return result | {"status": "run_full_mapping"}
        valid = f_b > np.maximum(g_b, e_b)
        from_g = (t_b[valid] - g_b[valid]) / (f_b[valid] - g_b[valid])
        from_e = (t_b[valid] - e_b[valid]) / (f_b[valid] - e_b[valid])
        low = np.quantile(np.minimum(from_g, from_e), 0.025)
        high = np.quantile(np.maximum(from_g, from_e), 0.975)
        point = (test - bg) / (f - bg)
    else:
        # Enforce sum(p)=1 while retaining signed estimates near the detection floor.
        # Populations satisfy d = M p, with rows as settings and columns as prep labels.
        def solve(data):
            design = data[..., :2] - data[..., 2:3]
            rhs = data[..., 3] - data[..., 2]
            p01 = (np.linalg.pinv(design) @ rhs[..., None])[..., 0]
            return 1 - p01.sum(axis=-1)

        design = mean[:, :2] - mean[:, 2:3]
        singular = np.linalg.svd(design, compute_uv=False)
        result["mapping_condition_number"] = float(singular[0] / singular[-1]) if singular[-1] else None
        if singular[-1] <= 0.05 or singular[0] / singular[-1] > 10:
            return result | {"status": "ill_conditioned_mapping"}
        point = solve(mean)
        low, high = np.quantile(solve(samples), [0.025, 0.975])

    status = "resolved" if low > 0 else "not_resolved"
    upper_end = float(min(high, 1))
    if high < 0 or low > 1:
        status = "inconsistent_with_population_model"
        upper_end = None
    return result | {
        "status": status,
        "p2": float(point),
        "interval_95": [float(low), float(high)],
        "upper_end_95_clipped": upper_end,
    }


def fetch_bright_fractions(handles, analyses, blocks, timeout):
    """Read IQCC 0.20 and native QOP result streams using their respective APIs."""
    from iqcc_cloud_client.qmm_cloud import CloudResultHandles

    is_cloud = isinstance(handles, CloudResultHandles)
    # IQCC execute returns the completed payload synchronously. Its wait method
    # is a no-op and its result wrapper does not accept flat_struct.
    if not is_cloud and not handles.wait_for_all_values(timeout=timeout):
        raise RuntimeError("Acquisition did not finish with all result streams.")
    bright = np.empty((blocks, len(analyses), len(STATES)))
    for a, analysis in enumerate(analyses):
        for s, preparation in enumerate(STATES):
            name = f"{analysis}_{preparation}"
            if name not in handles.keys():
                raise RuntimeError(f"Missing acquisition stream: {name}.")
            handle = handles.get(name)
            raw = handle.fetch_all() if is_cloud else handle.fetch_all(flat_struct=True)
            if isinstance(raw, dict) and "value" in raw:
                raw = raw["value"]
            values = np.asarray(raw)
            if values.dtype.names and "value" in values.dtype.names:
                values = values["value"]
            if values.size != blocks or (not is_cloud and handle.has_dataloss()):
                raise RuntimeError(f"Incomplete acquisition for {name}.")
            bright[:, a, s] = values.reshape(-1)
    return bright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qubit", required=True)
    parser.add_argument("--quam-state", type=Path, required=True)
    parser.add_argument("--test-op", default="x180")
    parser.add_argument("--f12-hz", type=float, required=True, help="Absolute EF RF frequency in Hz")
    parser.add_argument("--ref-amp", type=float, required=True, help="Calibrated GE reference cosine amplitude")
    parser.add_argument("--ref-ns", type=int, required=True)
    parser.add_argument("--ef-amp", type=float, required=True, help="Calibrated EF cosine amplitude")
    parser.add_argument("--ef-ns", type=int, required=True)
    parser.add_argument("--blocks", type=int, default=20)
    parser.add_argument("--shots-per-block", type=int, default=500)
    parser.add_argument("--timeout", type=int, default=180, help="Acquisition timeout in seconds")
    parser.add_argument("--background-tolerance", type=float, default=0.005, help="Allowed dark-fraction difference")
    parser.add_argument("--full", action="store_true", help="Acquire all three analyzer settings")
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--check", action="store_true", help="Build program and config locally")
    output.add_argument("--simulate", type=int, metavar="CLOCK_CYCLES", help="Request a QOP waveform simulation")
    output.add_argument("--run", type=Path, metavar="OUTPUT_JSON", help="Acquire hardware data into a new file")
    args = parser.parse_args()

    if (not np.isfinite(args.f12_hz) or args.f12_hz <= 0
            or any(not 0 < amp <= 0.25 for amp in (args.ref_amp, args.ef_amp))
            or any(ns < 16 or ns % 4 for ns in (args.ref_ns, args.ef_ns))
            or not 0 < args.background_tolerance < 1 or args.blocks < 8 or args.shots_per_block < 1
            or args.timeout < 1
            or (args.simulate is not None and args.simulate < 1)):
        parser.error("Invalid frequency, amplitude, duration, shot count, tolerance, or simulation duration.")
    if args.run and (args.run.exists() or not args.run.parent.is_dir()):
        parser.error("Choose a new output file in an existing directory.")

    from quam_config import Quam

    machine = Quam.load(args.quam_state)
    qubit = machine.qubits[args.qubit]
    experiment, config, analyses, ge_if, ef_if = build_program(machine, qubit, args)
    summary = {
        "qubit": args.qubit, "test_operation": args.test_op, "quam_state": str(args.quam_state),
        "f01_hz": qubit.xy.RF_frequency, "f12_hz": args.f12_hz, "ge_if_hz": ge_if, "ef_if_hz": ef_if,
        "reference_ge": {"amplitude": args.ref_amp, "duration_ns": args.ref_ns},
        "reference_ef": {"amplitude": args.ef_amp, "duration_ns": args.ef_ns},
        "thermalization_ns": qubit.thermalization_time, "state_order": STATES,
        "analysis_order": analyses, "blocks": args.blocks, "shots_per_block": args.shots_per_block,
        "background_tolerance": args.background_tolerance, "timeout_s": args.timeout,
    }
    if args.check:
        print(json.dumps(summary, indent=2))
        return

    qmm = machine.connect()
    if args.simulate is not None:
        from qm import SimulationConfig

        if not callable(getattr(qmm, "simulate", None)):
            parser.error("This manager does not expose QOP simulation. Use --check or a simulation-enabled QUAM state.")
        job = qmm.simulate(config, experiment, SimulationConfig(duration=args.simulate))
        job.wait_until("Done")
        print(job.get_simulated_waveform_report().to_string())
        return

    from iqcc_cloud_client.qmm_cloud import CloudQuantumMachinesManager
    from qualang_tools.multi_user import qm_session

    with qm_session(qmm, config) as qm:
        if isinstance(qmm, CloudQuantumMachinesManager):
            job = qm.execute(experiment, options={"timeout": args.timeout})
        else:
            job = qm.execute(experiment)
        bright = fetch_bright_fractions(job.result_handles, analyses, args.blocks, args.timeout)

    result = estimate(bright, args.background_tolerance)
    report = summary | {"bright_by_block": bright.tolist(), "estimate": result}
    with args.run.open("x") as output_file:
        json.dump(report, output_file, indent=2, allow_nan=False)
        output_file.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
