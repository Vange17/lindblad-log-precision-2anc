# Thesis numerical figures

This directory contains the figures used in the revised thesis and their
rebuilding script. The archived upstream data and the supplementary variance
results are both included in this repository; no download or new sampling is
needed to rebuild the plots.

- `gate.png`: thesis Figure 4.4, the 20-qubit gate-count comparison, unchanged
  by the variance update.
- `noise.png`: thesis Figure 4.5, the five-qubit observable-estimation errors
  with supplementary one-SEM ranges for the compensated protocol.
- `rebuild_numerical_figures.py`: plotting code for both figures.
- `test_rebuild_numerical_figures.py`: data alignment, SEM normalization,
  absolute-error mapping, and preservation-of-points tests.

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
