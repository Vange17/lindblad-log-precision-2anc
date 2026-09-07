# Thesis numerical figures

This directory contains the figures used in the revised thesis and their
rebuilding script. The archived upstream data and the supplementary variance
results are both included in this repository; no download or new sampling is
needed to rebuild the plots.

- `gate.png`: thesis Figure 4.4, the 20-qubit comparison of pre-synthesis CNOT
  counts and estimated T counts. The numerical values are unchanged; the
  rotation-related curves are now labeled T, not Rz.
- `noise.png`: thesis Figure 4.5, the five-qubit observable-estimation errors
  with supplementary one-SEM ranges for the compensated protocol.
- `rebuild_numerical_figures.py`: plotting code for both figures.
- `test_rebuild_numerical_figures.py`: gate-count units and layout, data alignment,
  SEM normalization, absolute-error mapping, and preservation-of-points tests.

## Figure 4.4 counting convention

The archived MATLAB arrays named `miti_Rz` and `trotter_Rz` already multiply
the rotation counts by 66. The thesis figure retains these values and reports
them as estimated T counts, using the convention in
[Zeng et al., PRX Quantum 6, 010359 (2025), Appendix G, Eq. (G1)](https://journals.aps.org/prxquantum/pdf/10.1103/PRXQuantum.6.010359#page=44):
the average T cost per Rz rotation is approximated by
`1.149 * log2(1 / delta) + 9.2`, giving approximately 66 at single-rotation
synthesis precision `delta = 1e-15`.
This estimate is based on the repeat-until-success synthesis method of
[Bocharov, Roetteler, and Svore](https://doi.org/10.1103/PhysRevLett.114.080502),
not a deterministic per-rotation gate count.

The per-step rotation counts before synthesis are 42 for the baseline and 51
for the compensated method. CNOT counts and the plotted target simulation
precision are also evaluated before rotation synthesis. Additional synthesis
CNOTs, ancillas, measurements, and accumulated synthesis errors are not included;
the figure is not a complete fault-tolerant resource analysis. The upstream
MATLAB source, its historical variable names and labels, and the archived MAT
file are preserved unchanged. Only the maintained plotting code uses the
corrected T-count names and labels. All plotted values, fits, and layout
settings are retained.

## Figure 4.5 supplementary statistics

Figure 4.5 retains the original 5,000,000-trajectory means. Its SEM estimates
come from 10,000 independent supplementary trajectories of the corrected
sampler, not from a new 5,000,000-trajectory run. The signed mean's one-SEM
range is mapped to the absolute-error axis, with endpoints `max(0, e-h)` and
`e+h`, where `e` is the plotted absolute error and `h` is the supplementary SEM.
These ranges are not confidence intervals or standard deviations of absolute
error. Applicability to the archived means assumes matching estimator
distributions; see the [variance supplement](../variance/README.md) for the
initialization correction and the unresolved historical provenance.

All 76 time points are available in the variance CSV. The figure displays the
same 39 selected points as the thesis plot before the SEM update. Figure sizes,
fonts, axes, and archived point values are preserved.

## Rebuild

From the repository root, with NumPy, SciPy, and Matplotlib available:

```sh
python3 figures/rebuild_numerical_figures.py miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots.mat --output-dir figures/generated
```

This rebuilds both figures into an ignored output directory, leaving the
committed images unchanged. Add `--figure error` to rebuild only Figure 4.5,
or `--figure gate` for Figure 4.4. A different variance table can be supplied
with `--variance-summary`; mismatched times, archived errors, or SEM
normalization are rejected.

The recorded exports use Times New Roman, with STIX math text. Install that
font separately if exact typography is required; Matplotlib may otherwise
substitute a locally available serif font. The recorded Python environment
uses NumPy 1.23.5, SciPy 1.13.1, and Matplotlib 3.8.2.

## Check

```sh
MPLBACKEND=Agg python3 -m unittest discover -s figures -p 'test_rebuild_numerical_figures.py' -v
```

Tests use the recorded benchmark data. They do not run new Monte Carlo
trajectories or replace the archived means.
