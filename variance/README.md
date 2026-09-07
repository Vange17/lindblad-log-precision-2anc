# Supplementary trajectory variance and standard errors

This directory provides a reproducible variance estimate for the five-qubit
real-time benchmark. It contains two independent batches of 5,000 trajectories,
their per-trajectory outputs at all 76 time points, and pooled statistics for
10,000 trajectories. The preliminary 200-trajectory pilot is not included.

The upstream MATLAB script and its original data file are unchanged. These
sampling and aggregation scripts do not regenerate figures; the separate
[figure rebuilding script](../figures/README.md) uses the recorded results.

## Statistics

For each time point, the sampled quantity is the exact observable expectation
along a random compensation trajectory, multiplied by its cumulative
normalization weight. It corresponds to `overlap_rhoNM2` in the upstream script.
There is no additional sampling of quantum measurement outcomes, so the
statistics do not include measurement-shot noise.

The sample variance uses the pooled 10,000 outputs and denominator 9,999. The
trajectory standard deviation is its square root. The standard error of a mean
(SEM) is the trajectory standard deviation divided by the square root of the
number of independent trajectories used in that mean.

| Time | Trajectory variance | Trajectory standard deviation | SEM extrapolated to 5,000,000 trajectories |
| --- | --- | --- | --- |
| 0.5 | 0.34476 | 0.58716 | 0.00026259 |
| 1.0 | 0.38995 | 0.62446 | 0.00027927 |
| 1.5 | 0.64699 | 0.80435 | 0.00035972 |

The last column is an extrapolation from the supplementary variance estimate,
not a new run of 5,000,000 trajectories. These are standard errors, not confidence
intervals. No bootstrap calculation or confidence multiplier is applied.

All time points share the same set of independent trajectories; different time
points within one trajectory are not independent samples. The reported SEM
belongs to the signed observable estimator, not to its absolute deviation from
the reference. In particular, it is not the standard deviation of absolute error.

In thesis Figure 4.5, the historical mean estimates are retained. The supplementary
one-SEM range of the signed deviation is mapped to the absolute-error axis:
for an absolute error `e` and SEM `h`, the displayed endpoints are
`max(0, e-h)` and `e+h`. These ranges are supplementary estimates subject to the
estimator-distribution limitation below, not confidence intervals or empirical
error bars recovered from the historical trajectories.

The variance estimate itself has sampling uncertainty. An approximate standard
error based on the sample fourth central moment is retained as a diagnostic.
Its relative values at the three times above are approximately 1.15%, 2.74%,
and 5.49%; the last digits of the variance estimates should not be treated as
high-precision results.

## Data and provenance

- Upstream: [xiaogangli1169/lindblad-log-precision-2anc](https://github.com/xiaogangli1169/lindblad-log-precision-2anc).
- Upstream revision: `db4e49f9ab20aaf907803dc788f71206efbdfd6a`.
- Input: the original `miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots.mat` in the repository root.
- Input MAT Git blob SHA-1: `29c0648dfe081f414d15f0fc864bb7d13223782c`.
- Model: five qubits, time step 0.02, times 0 through 1.5, truncation order 7. Hamiltonian, dilation, initial state, observable, compensation coefficients, and sampling probabilities are loaded directly from the original MAT file.
- Random streams: MATLAB `twister`, seeds 20260908 and 20260909, respectively.
- Recorded environment: MATLAB R2023b Update 2; aggregation with Python 3.9, NumPy 1.23.5 and SciPy 1.13.1.
- Sampling performed on 2026-09-07, using a single MATLAB computation thread on an Apple M2 with 16 GiB RAM. The two batches took approximately 55 and 57 seconds, excluding startup.

The new sampler initializes `xlab_prod=1` before the `kA` branch at every time
step. The published upstream script initializes it only inside the `kA==0`
branch, leaving it undefined on some first-step draws or carrying it from an
earlier step. This correction is included only in `run_variance.m`.

Consequently, the supplementary data describe the corrected sampler. They do
not recover the empirical variance of the historical samples. Whether the
archived upstream means were generated with the exact published initialization
logic is not established here. Applying the extrapolated SEM to those means
requires matching the estimator distributions.

## Files

| File | Contents |
| --- | --- |
| `run_variance.m` | Standalone sequential sampler; repository-relative input path; no extra MATLAB toolboxes required |
| `summarize_variance.py` | Validates and pools the two batches, then exports variance, standard deviation and SEM |
| `verify_variance.m` | Small initialization, numerical equivalence and overwrite-protection checks |
| `results/batch_a_5000.mat`, `results/batch_b_5000.mat` | Recorded independent batches; `values` is a 5,000-by-76 matrix of reweighted outputs |
| `results/variance_summary.csv` | Statistics at all 76 time points |
| `results/variance_summary.json` | Counts, seeds, software versions, hashes, upstream revision and recorded checks |
| `tests/reference_20.mat` | Twenty trajectories computed using the original matrix-exponential and ancilla-dilation operations, with the initialization correction |
| `tests/recorded_validation.mat` | Initial regression and equivalence results recorded with the supplementary runs |

The main plotting quantity for a 5,000,000-trajectory mean is the CSV column
`mean_se_extrapolated_N5000000`. It must not be confused with
`trajectory_std_ddof1`, `mean_se_N10000`, or `trajectory_variance_se_plugin`.

## Reproduce the pooled statistics

From the repository root, with NumPy and SciPy already available:

```sh
python3 variance/summarize_variance.py
```

This reads the recorded batches and writes new summary files to
`variance/generated/`, leaving `results/` unchanged. Reusing an existing output
file is rejected; choose a new directory with `--output-dir` when needed.

## Reproduce the trajectories

From the repository root in MATLAB:

```matlab
addpath('variance');
mkdir(fullfile('variance','generated','rerun'));
run_variance(5000,20260908,fullfile('variance','generated','rerun','batch_a_5000'));
run_variance(5000,20260909,fullfile('variance','generated','rerun','batch_b_5000'));
```

Then aggregate the newly generated batches:

```sh
python3 variance/summarize_variance.py --data-dir variance/generated/rerun --output-dir variance/generated/rerun-summary
```

The sampler accepts 2 through 10,000 trajectories per call and refuses to
overwrite existing output files. Each MAT output retains the reweighted
per-trajectory values. A local inverse-CDF categorical sampler replaces the
upstream `randsrc`, preserving the sampling probabilities without requiring the
Communications Toolbox. It does not reproduce the historical random sequence.

The Hamiltonian step is precomputed, and the ancilla partial trace is evaluated
with two equivalent Kraus terms. The recorded maximum difference against the
reference matrix operations is 4.44e-16 per trajectory; the deterministic
baseline agrees with the original archive to 1.11e-15. The original sampler's
large symbolic preprocessing and multi-gigabyte preallocations are not needed.

## Quick checks

The following runs only 22 trajectories, not the complete supplementary experiment:

```matlab
addpath('variance');
mkdir(fullfile('variance','generated','quick-check'));
verify_variance(fullfile('variance','generated','quick-check'));
```

The seed-1989 test exercises the first-step branch that requires the initialization
correction. The seed-20260907 test compares every output of 20 trajectories with
the recorded reference. Variance/SEM identities, the deterministic baseline and
refusal to overwrite outputs are also checked.
