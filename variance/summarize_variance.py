"""Pool two fixed-size independent runs; keep the 200-trajectory pilot separate."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
REVISION = 'db4e49f9ab20aaf907803dc788f71206efbdfd6a'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=ROOT / 'results')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'generated')
    args = parser.parse_args()
    output_dir = args.output_dir
    for name in ('variance_summary.csv', 'variance_summary.json'):
        if (output_dir / name).exists():
            raise FileExistsError(f'Output exists; choose a new output directory: {output_dir / name}')
    source_hash = hashlib.sha256((REPO / 'Miti_Ising_OQS4_noshots_github.m').read_bytes()).hexdigest()
    if source_hash != '44dba804e95d86ea9c8262b197240740cba63dc641cbfa7c34dc6e5a0ac31066':
        raise ValueError('Upstream MATLAB source does not match the pinned revision.')
    batch_paths = [args.data_dir / f'batch_{label}_5000.mat' for label in ('a', 'b')]
    batch_hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in batch_paths}
    batches = [loadmat(path, squeeze_me=True) for path in batch_paths]
    validation = loadmat(ROOT / 'tests' / 'recorded_validation.mat', squeeze_me=True)
    assert bool(validation['caught'])
    assert float(validation['trajectory_max_difference']) < 1e-11
    for batch, seed in zip(batches, (20260908, 20260909)):
        assert int(batch['num_NMt']) == 5000 and int(batch['seed']) == seed
        assert batch['values'].shape == (5000, 76)
        assert np.isfinite(batch['values']).all() and np.isrealobj(batch['values'])
        assert float(batch['baseline_max_difference']) < 1e-11
        np.testing.assert_allclose(np.var(batch['values'], axis=0, ddof=1),
                                   batch['trajectory_variance'], rtol=1e-12, atol=1e-14)
    for key in ('t_vec', 'mu_step', 'exact_mean', 'archive_mean'):
        np.testing.assert_array_equal(batches[0][key], batches[1][key])
    times = batches[0]['t_vec']
    np.testing.assert_allclose(times, np.arange(76) * 0.02, atol=1e-14)
    values = np.concatenate([b['values'] for b in batches])
    n = len(values)
    assert np.all(values[:, 0] == 1)
    assert np.all(np.abs(values) <= float(batches[0]['mu_step'])**np.arange(76) + 1e-10)
    means = values.mean(axis=0)
    variances = values.var(axis=0, ddof=1)
    pooled_m2 = sum(4999 * b['trajectory_variance']
                    + 5000 * (b['trajectory_mean'] - means)**2 for b in batches)
    np.testing.assert_allclose(variances, pooled_m2 / (n - 1), rtol=1e-12, atol=1e-14)
    fourth = np.mean((values - means)**4, axis=0)
    # Plug-in standard error of unbiased sample variance, not an exact CI.
    variance_of_variance = (fourth - (n - 3) / (n - 1) * variances**2) / n
    assert variance_of_variance.min() >= -1e-14
    variance_se = np.sqrt(np.maximum(variance_of_variance, 0))
    se_run = np.sqrt(variances / n)
    se_original_n = np.sqrt(variances / 5_000_000)
    z = np.divide(means - batches[0]['exact_mean'], se_run,
                  out=np.zeros_like(means), where=se_run > 0)
    relative_variance_se = np.divide(variance_se, variances,
                                     out=np.zeros_like(means), where=variances > 0)
    columns = {
        'time': times,
        'trajectory_mean_N10000': means,
        'trajectory_variance_ddof1': variances,
        'trajectory_std_ddof1': np.sqrt(variances),
        'trajectory_variance_se_plugin': variance_se,
        'mean_se_N10000': se_run,
        'mean_se_extrapolated_N5000000': se_original_n,
        'exact_mean': batches[0]['exact_mean'],
        'archive_mean_N5000000': batches[0]['archive_mean'],
        'variance_batch_a_N5000': batches[0]['trajectory_variance'],
        'variance_batch_b_N5000': batches[1]['trajectory_variance'],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / 'variance_summary.csv').open('x', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))
    summary = {
        'source_url': 'https://github.com/xiaogangli1169/lindblad-log-precision-2anc',
        'source_revision': REVISION,
        'source_script_sha256': source_hash,
        'input_batch_sha256': batch_hashes,
        'reproduction_script_sha256': hashlib.sha256((ROOT / 'run_variance.m').read_bytes()).hexdigest(),
        'analysis_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_mat_git_blob_sha1': '29c0648dfe081f414d15f0fc864bb7d13223782c',
        'statistic': 'Unbiased trajectory variance (ddof=1); trajectory SD; SEM at the specified sample count. No confidence intervals.',
        'supplementary_trajectory_count': n,
        'matlab_version': str(batches[0]['software_version']),
        'numpy_version': np.__version__,
        'scipy_version': scipy.__version__,
        'seeds': [20260908, 20260909],
        'rng': 'MATLAB twister; inverse-CDF categorical sampling',
        'trajectories_per_batch': 5000,
        'pooled_trajectories': n,
        'pilot_included': False,
        'original_trajectory_count_for_se_extrapolation': 5_000_000,
        'mu_step': float(batches[0]['mu_step']),
        'calculation_seconds': [float(b['elapsed_seconds']) for b in batches],
        'initialization_regression_seed': int(validation['rare_seed']),
        'reference_trajectory_max_difference': float(validation['trajectory_max_difference']),
        'baseline_max_difference': float(batches[0]['baseline_max_difference']),
        'max_absolute_z_to_exact_mean': float(np.max(np.abs(z))),
        'max_relative_variance_se_plugin_nonzero_times': float(relative_variance_se[1:].max()),
        'mean_se_N5000000_min_nonzero_time': float(se_original_n[1:].min()),
        'mean_se_N5000000_max': float(se_original_n.max()),
        'interpretation': 'Independent supplementary variance estimate with per-step xlab_prod reset; not the empirical variance of the historical samples. SEM for N=5000000 is extrapolated, not a new N=5000000 run. Exact trace per trajectory; no measurement-shot noise.',
        'selected_times': [
            dict(time=float(times[j]), variance=float(variances[j]),
                 trajectory_std=float(np.sqrt(variances[j])),
                 variance_se_plugin=float(variance_se[j]),
                 mean_se_N5000000=float(se_original_n[j]),
                 variance_batch_a=float(batches[0]['trajectory_variance'][j]),
                 variance_batch_b=float(batches[1]['trajectory_variance'][j]),
                 mean_z_to_exact=float(z[j]))
            for j in (0, 25, 50, 75)
        ],
    }
    with (output_dir / 'variance_summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
