# Measuring Leakage After a Single X01 Pulse on Qolab

Estimate the population in the second excited state, **P₂**, immediately after one calibrated X₀₁(π) pulse applied to a thermally reset transmon. Start with Q3 and repeat the same experiment on each available qubit. The measurement samples final population after one pulse.

## Calibrate the reference pulses

1. At the experimental flux bias, confirm f₀₁ with `03a_qubit_spectroscopy` and the target X₀₁(π) operation with `04b_power_rabi`. Confirm the standard binary readout with `07_iq_blobs`.
2. Calibrate a separate, spectrally selective X₀₁ reference π pulse. The script constructs it as a plain cosine pulse at the configured f₀₁. Choose a longer duration than the gate under test and determine its amplitude by a power Rabi scan with that exact envelope. Check its preparation and transfer contrast independently. Using the test pulse as the reference could absorb the leakage of interest into the calibration.
3. Find f₁₂ with the sequence reset → reference X₀₁ → swept EF drive → reference X₀₁ → readout. A near-resonant EF drive transfers population into |2⟩, which is largely unaffected by the final GE pulse. Compare against matched no-drive and off-resonant controls. At the observed resonance, sweep the EF cosine amplitude with the same sequence and identify the first transfer maximum and the return near twice its pulse area. Use the measured frequency and π amplitude at the selected pulse duration. For Q3, prior spectroscopy puts the resonance near f₀₁ − 303 MHz. Find the resonance separately for each remaining qubit.

The Python runner takes those independently calibrated amplitudes and durations as inputs. Acquire calibration scans using the repository's spectroscopy and power-Rabi sweep patterns. The EF pulse uses the existing XY element at a temporary intermediate frequency. Confirm that the required IF and waveform lie within the backend's supported range. The runner checks local program construction, while QOP simulation checks controller scheduling and waveforms.

## Acquire four interleaved preparations

Each shot uses a thermal reset. The test branch then plays the target X₀₁ exactly once. All four branches subsequently apply the same analyzer, **X₁₂(π) → reference X₀₁(π) → binary readout**.

| Branch | Preparation after reset | Ideal state after the analyzer |
|---|---|---|
| `g` | No pulse | \(\lvert 1\rangle\) |
| `e` | Reference X₀₁ | \(\lvert 2\rangle\) |
| `f` | Reference X₀₁ → X₁₂ | \(\lvert 0\rangle\) |
| `test` | One target X₀₁ | Depends on the test pulse |

The analyzer maps pre-readout |2⟩ onto |0⟩. Analyze immediately after preparation and keep the pulse timing and readout the same across branches. The |2⟩ preparation is the positive control. The |0⟩ and |1⟩ preparations measure how often the analyzer and readout produce a dark result in the absence of intended leakage.

Start with 20 acquisition blocks of 500 shots per branch. Every block contains all four preparations. Keep blocks separate for analysis, and repeat the run if the contrast or the uncertainty changes substantially over time.

## Estimate and report P₂

If the two background dark fractions, d₀ and d₁, agree within the required tolerance and the |2⟩ control has clearly higher dark contrast, take d_bg = (d₀ + d₁)/2 and estimate

```text
P2 = (d_test - d_bg) / (d_2 - d_bg)
```

The runner checks background equivalence using a block bootstrap. The default allowable background difference is 0.005 in absolute dark fraction. Set a tighter tolerance when a 0.5 percentage-point difference divided by the measured contrast would be too large for the experiment. The reported interval combines block-bootstrap variation with the range obtained from the two background choices.

If the backgrounds differ or equivalence has insufficient precision, rerun with `--full`. This adds no analyzer and reference X₀₁-only settings to the same acquisition. For each setting the three reference preparations give one row of the response matrix M. The runner fits d = M p with p₀ + p₁ + p₂ = 1 and checks the conditioning. The fit retains signed estimates near zero to expose baseline bias. If the positive control gives unresolved contrast, check f₁₂, both reference pulses and the readout before interpreting the test branch.

Report the qubit, test operation, f₀₁, f₁₂, shots per branch, all four dark fractions, P₂ and its interval. For an unresolved result, quote the upper endpoint of the nominal 95% interval and the reference quality. The interval reflects acquisition and calibration counting variation across blocks. Residual preparation errors, relaxation during mapping, and analyzer sensitivity to coherences can add systematic bias. Where those effects could match the claimed leakage, report the result as conditional on the calibration and verify it with an independent reference or a phase-cycled analyzer. Population above |2⟩ requires additional controls.

## Run from the Qolab Start checkout

Use the existing `qolab-start` virtual environment. Set the shell variables below from the calibrations for the selected qubit. The state path must point to the active QUAM directory containing `state.json` and `wiring.json`.

```bash
STATE=qualibration_graphs/superconducting/calibrations/1Q_calibrations/quam_state
SRC=../project/src/single_x01_leakage.py

.venv/bin/python "$SRC" --quam-state "$STATE" --qubit Q3 \
  --f12-hz "$F12_HZ" --ref-amp "$REF_AMP" --ref-ns "$REF_NS" \
  --ef-amp "$EF_AMP" --ef-ns "$EF_NS" --check
```

The active Cloud QUAM state selects the IQCC Cloud Quantum Machines Manager. Configure the account token using `iqcc-cloud setup` and confirm access to the assigned backend before collecting data. QUAlibrate can still be used for prerequisite calibrations and state management. The standalone runner submits its QUA program through `Quam.connect()`.

To collect hardware data after verifying the settings, use `--run Q3_single_x01.json`. The JSON stores block-wise bright fractions, run settings and the estimate. Use a new filename per qubit and run. The full-mapping version collects three times as many shots as the default at the same number of blocks. Change `--timeout` when allowing more time for a long acquisition.

Use `--simulate CLOCK_CYCLES` with a simulation-enabled QOP manager to inspect the waveform report. Use enough cycles to cover the thermal wait and the branches of interest. The IQCC Cloud 0.20.0 manager in this checkout does not expose the QOP simulation method.
