"""Rebuild Figures 4.4 and 4.5 from the archived numerical data.

Data and the original MATLAB implementation:
https://github.com/xiaogangli1169/lindblad-log-precision-2anc
Pinned source revision: db4e49f9ab20aaf907803dc788f71206efbdfd6a

The MATLAB file stores several symbolic arrays that scipy cannot decode.  The
gate counts are therefore recomputed from the deterministic formulas in
Miti_Ising_OQS4_noshots_github.m.  Its factor of 66 is retained as an estimated
T-gate cost per R_Z rotation, following PRX Quantum 6, 010359, Eq. (G1), at
single-rotation synthesis precision 1e-15. CNOT counts remain pre-synthesis;
additional synthesis overheads and accumulated synthesis errors are not included.
The real-time errors are read directly from
the archived numerical arrays.  Figure 4.5 retains these means and adds SEM
estimates from independent supplementary trajectories of the corrected sampler;
matching its estimator distribution to the historical run is not established.
"""

import argparse
import io
import math
import struct
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, whosmat
from scipy.optimize import brentq


SOURCE_REVISION = "db4e49f9ab20aaf907803dc788f71206efbdfd6a"
SOURCE_URL = "https://github.com/xiaogangli1169/lindblad-log-precision-2anc"
VARIANCE_SUMMARY = (
    Path(__file__).resolve().parents[1]
    / "variance/results/variance_summary.csv"
)
EXPECTED_VARIABLES = {
    "A_abs",
    "deltaJ0_noshots",
    "deltaNM_noshots",
    "epslog10",
    "t_vec",
}

BLUE = "#0072B2"
ORANGE = "#D55E00"


def load_numeric_variables(mat_path, names):
    """Read selected numeric arrays from a MATLAB v5 file.

    The source archive contains MATLAB symbolic objects that make a whole-file
    scipy load fail.  MATLAB v5 stores each workspace variable in an independent
    top-level data element, so this function decodes only the requested numeric
    elements and deliberately skips symbolic ones.
    """

    data = mat_path.read_bytes()
    if not data.startswith(b"MATLAB 5.0 MAT-file"):
        raise ValueError(f"Expected a MATLAB v5 file: {mat_path}")

    header = data[:128]
    offset = 128
    values = {}

    while offset + 8 <= len(data) and set(values) != set(names):
        data_type, byte_count = struct.unpack_from("<II", data, offset)
        # MATLAB does not pad top-level miCOMPRESSED elements, although other
        # data element types are aligned to 64-bit boundaries.
        stored_count = byte_count if data_type == 15 else (byte_count + 7) // 8 * 8
        end = offset + 8 + stored_count
        if end > len(data):
            raise ValueError(f"Truncated MATLAB data element at byte {offset}")

        element = data[offset:end]
        try:
            variables = whosmat(io.BytesIO(header + element))
        except (OSError, TypeError, ValueError):
            variables = []

        if variables and variables[0][0] in names:
            variable_name = variables[0][0]
            loaded = loadmat(io.BytesIO(header + element))
            values[variable_name] = np.asarray(loaded[variable_name]).squeeze()

        offset = end

    missing = set(names) - set(values)
    if missing:
        raise ValueError(f"Missing MATLAB variables: {sorted(missing)}")
    return values


def compute_gate_counts(variables):
    """Return pre-synthesis CNOT counts and estimated T counts for 20 qubits."""

    a2 = float(variables["A_abs"][1])
    gamma = 1.5
    dissipator_norm = math.sqrt(gamma)
    fitted_lambda = 1.7
    total_time = 1.0
    number_of_spins = 20
    hamiltonian_norm = 1.5 * number_of_spins
    lindbladian_norm = hamiltonian_norm + dissipator_norm**2

    def normalization_log(time_step):
        return total_time * (
            2 * 3.5 * (2 * hamiltonian_norm) ** 4 * time_step**3
            + (
                2 * a2 * dissipator_norm**4
                + 16 * lindbladian_norm**2
            )
            * time_step
        )

    compensated_time_step = brentq(
        lambda time_step: normalization_log(time_step) - 3,
        1e-8,
        1,
    )
    compensated_steps = math.ceil(total_time / compensated_time_step)

    def compensated_error(order):
        hamiltonian_base = (
            2 * math.e * hamiltonian_norm * compensated_time_step / (order + 1)
        )
        hamiltonian_error = (
            2 * hamiltonian_base ** (order + 1)
            + hamiltonian_base ** (2 * order + 2)
        )
        dissipative_error = (
            2
            * a2
            * fitted_lambda ** (order - 1)
            * dissipator_norm ** (2 * order + 2)
            * compensated_time_step ** (order + 1)
        )
        trotter_error = (
            4
            * math.e
            * lindbladian_norm
            * compensated_time_step
            / (order + 1)
        ) ** (order + 1)
        return hamiltonian_error + dissipative_error + trotter_error

    eps_log10 = np.asarray(variables["epslog10"], dtype=float)
    precision = 10**eps_log10
    truncation_orders = np.array(
        [
            brentq(
                lambda order: compensated_error(order) - target_precision,
                0,
                100,
            )
            for target_precision in precision
        ]
    )

    baseline_coefficient = (
        0.04 * (number_of_spins - 1)
        + 2 * a2 * gamma**2
        + (4 * math.e * (0.3 * number_of_spins + gamma) / 2) ** 2
    )
    baseline_time_step = precision / baseline_coefficient
    baseline_steps = np.ceil(total_time / baseline_time_step)

    compensated_cnot = (
        2 * number_of_spins + 8 + 16 * np.ceil(truncation_orders)
    ) * compensated_steps
    compensated_cnot_continuous = (
        2 * number_of_spins + 8 + 16 * truncation_orders
    ) * compensated_steps
    # 66 is the average T cost per rotation, not a count of R_Z gates.
    # Zeng et al., https://doi.org/10.1103/PRXQuantum.6.010359, Eq. (G1).
    compensated_t = np.full_like(
        precision,
        66 * (2 * number_of_spins + 11) * compensated_steps,
    )
    baseline_cnot = (2 * number_of_spins + 4) * baseline_steps
    baseline_t = 66 * (2 * number_of_spins + 2) * baseline_steps

    return {
        "eps_log10": eps_log10,
        "precision": precision,
        "baseline_t": baseline_t,
        "baseline_cnot": baseline_cnot,
        "compensated_t": compensated_t,
        "compensated_cnot": compensated_cnot,
        "compensated_cnot_continuous": compensated_cnot_continuous,
    }


def configure_plot_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman"],
            "mathtext.fontset": "stix",
            "font.size": 20,
            "axes.labelsize": 20,
            "axes.linewidth": 2.5,
            "legend.fontsize": 18,
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )


def draw_gate_figure(gate_data, output_path):
    precision = gate_data["precision"]
    eps_log10 = gate_data["eps_log10"]
    fit_log_precision = np.linspace(eps_log10.min(), eps_log10.max(), 400)
    fit_precision = 10**fit_log_precision

    figure, axes = plt.subplots(figsize=(12.0, 5.85))
    figure.subplots_adjust(left=0.12, right=0.76, bottom=0.22, top=0.95)

    marker_options = {"linestyle": "none", "markersize": 10, "markeredgewidth": 1.0}
    axes.loglog(
        precision,
        gate_data["baseline_t"],
        marker="D",
        color=BLUE,
        label=r"$T$ (Trotter-like)",
        **marker_options,
    )
    axes.loglog(
        precision,
        gate_data["baseline_cnot"],
        marker="^",
        color=BLUE,
        label="CNOT (Trotter-like)",
        **marker_options,
    )
    axes.loglog(
        precision,
        gate_data["compensated_t"],
        marker="D",
        color=ORANGE,
        label=r"$T$ (Ours)",
        **marker_options,
    )
    axes.loglog(
        precision,
        gate_data["compensated_cnot"],
        marker="^",
        color=ORANGE,
        label="CNOT (Ours)",
        **marker_options,
    )

    fit_specs = (
        ("baseline_t", True),
        ("baseline_cnot", True),
        ("compensated_t", True),
        ("compensated_cnot_continuous", False),
    )
    for index, (key, log_fit) in enumerate(fit_specs):
        values = gate_data[key]
        if log_fit:
            coefficients = np.polyfit(eps_log10, np.log10(values), 1)
            fitted = 10 ** np.polyval(coefficients, fit_log_precision)
        else:
            coefficients = np.polyfit(eps_log10, values, 1)
            fitted = np.polyval(coefficients, fit_log_precision)
        axes.loglog(
            fit_precision,
            fitted,
            color="black",
            linestyle=(0, (4, 3)),
            linewidth=2.0,
            label="Fit" if index == 0 else None,
            zorder=1,
        )

    annotations = (
        (r"$\mathcal{O}(\varepsilon^{-1})$", 0.92),
        (r"$\mathcal{O}(\varepsilon^{-1})$", 0.79),
        (r"$\mathcal{O}(\varepsilon^{0})$", 0.20),
        (
            r"$\mathcal{O}(\log(1/\varepsilon))$",
            0.08,
        ),
    )
    for label, axes_y in annotations:
        axes.text(
            1.025,
            axes_y,
            label,
            transform=axes.transAxes,
            ha="left",
            va="center",
            fontsize=20,
            clip_on=False,
        )

    axes.set_xlim(1.25e-2, 7.5e-9)
    axes.set_ylim(1e5, 1e15)
    axes.set_xticks([1e-2, 1e-4, 1e-6, 1e-8])
    axes.set_yticks([1e5, 1e7, 1e9, 1e11, 1e13, 1e15])
    axes.set_xlabel(r"Precision, $\varepsilon$")
    axes.set_ylabel("Gate count per sampled circuit")
    axes.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black")

    figure.savefig(
        output_path,
        dpi=330,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.02,
        metadata={
            "Source": f"{SOURCE_URL}/commit/{SOURCE_REVISION}",
            "GateCountConvention": (
                "Pre-synthesis CNOT counts; estimated T counts use 66 T gates per "
                "R_Z at synthesis precision 1e-15 (PRX Quantum 6, 010359, Eq. G1). "
                "Additional synthesis overheads and accumulated synthesis errors "
                "are not included."
            ),
        },
    )
    plt.close(figure)


def load_supplementary_sem(summary_path, time, compensated_error):
    """Load the corrected-sampler SEM, extrapolated to the archived sample count."""

    summary = np.genfromtxt(summary_path, delimiter=",", names=True, dtype=float)
    required = {
        "time",
        "trajectory_variance_ddof1",
        "mean_se_extrapolated_N5000000",
        "archive_mean_N5000000",
        "exact_mean",
    }
    if summary.shape != (76,) or not required.issubset(summary.dtype.names or ()):
        raise ValueError("Unexpected supplementary variance table shape or columns")
    if not all(np.isfinite(summary[name]).all() for name in required):
        raise ValueError("Non-finite supplementary variance data")
    if not np.allclose(summary["time"], time, rtol=0, atol=1e-12):
        raise ValueError("Supplementary times do not match the archived time grid")
    archive_error = np.abs(summary["archive_mean_N5000000"] - summary["exact_mean"])
    if not np.allclose(archive_error, compensated_error, rtol=1e-10, atol=1e-12):
        raise ValueError("Supplementary table does not match the archived errors")
    variance = summary["trajectory_variance_ddof1"]
    sem = summary["mean_se_extrapolated_N5000000"]
    if np.any(variance < 0) or np.any(sem < 0):
        raise ValueError("Negative supplementary variance or SEM")
    if not np.allclose(sem, np.sqrt(variance / 5_000_000), rtol=1e-12, atol=1e-14):
        raise ValueError("Supplementary SEM does not use the archived sample count")
    return sem


def draw_error_figure(variables, output_path, variance_summary=VARIANCE_SUMMARY):
    time = np.asarray(variables["t_vec"], dtype=float)
    baseline_error = np.asarray(variables["deltaJ0_noshots"], dtype=float)
    compensated_error = np.asarray(variables["deltaNM_noshots"], dtype=float)
    if not (time.shape == baseline_error.shape == compensated_error.shape == (76,)):
        raise ValueError("Unexpected real-time dataset shape")
    if not all(
        np.isfinite(data).all() for data in (time, baseline_error, compensated_error)
    ):
        raise ValueError("Non-finite real-time data")
    if np.any(baseline_error < 0) or np.any(compensated_error < 0):
        raise ValueError("Negative archived absolute error")
    sem = load_supplementary_sem(variance_summary, time, compensated_error)

    selected = np.r_[0:76:2, 75]
    figure, axes = plt.subplots(figsize=(8.0, 6.0))
    figure.subplots_adjust(left=0.18, right=0.97, bottom=0.20, top=0.91)

    axes.plot(
        time[selected],
        baseline_error[selected],
        marker="o",
        linestyle="none",
        markersize=10,
        color=BLUE,
        label="Trotter-like",
    )
    # Map the signed mean's +/- one-SEM interval to the absolute-error axis.
    # These are supplementary uncertainty scales, not SDs of absolute error.
    axes.errorbar(
        time[selected],
        compensated_error[selected],
        yerr=np.vstack(
            (np.minimum(compensated_error[selected], sem[selected]), sem[selected])
        ),
        marker="o",
        linestyle="none",
        markersize=10,
        color=ORANGE,
        elinewidth=1.0,
        capsize=2.5,
        capthick=1.0,
        label="Ours",
    )

    axes.set_xlim(0, 1.52)
    axes.set_ylim(0, 2.0e-3)
    axes.set_xticks([0, 0.5, 1.0, 1.5])
    axes.set_yticks(np.linspace(0, 2.0e-3, 5))
    axes.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3), useMathText=True)
    axes.set_xlabel(r"Simulation time, $T$")
    axes.set_ylabel("Observable estimation error")
    axes.legend(loc="upper right", frameon=True, fancybox=False, edgecolor="black")

    figure.savefig(
        output_path,
        dpi=330,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.02,
        metadata={
            "Source": f"{SOURCE_URL}/commit/{SOURCE_REVISION}",
            "SupplementarySEM": (
                "sqrt(s_T^2/5000000), with s_T^2 estimated from 10000 corrected-sampler "
                "trajectories. Archived means are unchanged; matching estimator "
                "distributions is unverified. One-SEM ranges are mapped to absolute error."
            ),
        },
    )
    plt.close(figure)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mat_file",
        type=Path,
        help="Path to miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots.mat",
    )
    parser.add_argument(
        "--figure",
        choices=("gate", "error", "both"),
        default="both",
        help="Select the figure(s) to rebuild; use error to leave Figure 4.4 untouched",
    )
    parser.add_argument(
        "--variance-summary",
        type=Path,
        default=VARIANCE_SUMMARY,
        help="Supplementary variance CSV used for Figure 4.5",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for gate.png and noise.png",
    )
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if not arguments.mat_file.is_file():
        raise FileNotFoundError(arguments.mat_file)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    configure_plot_style()
    variables = load_numeric_variables(arguments.mat_file, EXPECTED_VARIABLES)
    if arguments.figure in ("gate", "both"):
        gate_data = compute_gate_counts(variables)
        draw_gate_figure(gate_data, arguments.output_dir / "gate.png")
    if arguments.figure in ("error", "both"):
        draw_error_figure(
            variables, arguments.output_dir / "noise.png", arguments.variance_summary
        )


if __name__ == "__main__":
    main()
