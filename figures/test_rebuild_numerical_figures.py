"""Gate-count and supplementary SEM checks using the recorded benchmark data."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

import rebuild_numerical_figures as figures


class NumericalFigureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mat_path = figures.VARIANCE_SUMMARY.parents[2] / (
            "miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots.mat"
        )
        cls.variables = figures.load_numeric_variables(mat_path, figures.EXPECTED_VARIABLES)
        cls.summary = np.genfromtxt(
            figures.VARIANCE_SUMMARY, delimiter=",", names=True, dtype=float
        )
        figures.configure_plot_style()

    def test_gate_counts_use_estimated_t_units(self):
        counts = figures.compute_gate_counts(self.variables)
        self.assertNotIn("baseline_rz", counts)
        self.assertNotIn("compensated_rz", counts)
        np.testing.assert_array_equal(
            counts["baseline_t"], 66 * 42 * (counts["baseline_cnot"] / 44)
        )
        np.testing.assert_array_equal(counts["compensated_t"], np.full(13, 17826336))
        np.testing.assert_array_equal(
            counts["compensated_cnot"],
            [338944, 338944, 423680, 423680, 423680, 423680,
             508416, 508416, 508416, 508416, 508416, 593152, 593152],
        )
        self.assertEqual(counts["baseline_t"][0], 463043196)
        self.assertEqual(counts["baseline_t"][-1], 463041150890472)

    def test_gate_plot_labels_units_and_preserves_layout(self):
        counts = figures.compute_gate_counts(self.variables)
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(figures.plt, "close"):
                figures.draw_gate_figure(counts, Path(directory) / "gate.png")
                figure = figures.plt.gcf()
            try:
                axes = figure.axes[0]
                labels = [text.get_text() for text in axes.get_legend().get_texts()]
                self.assertEqual(labels, [
                    r"$T$ (Trotter-like)", "CNOT (Trotter-like)",
                    r"$T$ (Ours)", "CNOT (Ours)", "Fit",
                ])
                for line, key in zip(axes.lines[:4], (
                    "baseline_t", "baseline_cnot", "compensated_t", "compensated_cnot"
                )):
                    np.testing.assert_array_equal(line.get_xdata(), counts["precision"])
                    np.testing.assert_array_equal(line.get_ydata(), counts[key])
                np.testing.assert_array_equal(figure.get_size_inches(), [12, 5.85])
                np.testing.assert_allclose(axes.get_position().bounds, [0.12, 0.22, 0.64, 0.73])
                self.assertEqual(axes.get_xlim(), (1.25e-2, 7.5e-9))
                self.assertEqual(axes.get_ylim(), (1e5, 1e15))
                self.assertEqual(axes.get_ylabel(), "Gate count per sampled circuit")
                self.assertEqual(axes.yaxis.label.get_fontsize(), 20)
                self.assertEqual(len(axes.texts), 4)
                self.assertTrue(all(text.get_position()[0] == 1.025 for text in axes.texts))
            finally:
                figures.plt.close(figure)

    def load_sem(self):
        return figures.load_supplementary_sem(
            figures.VARIANCE_SUMMARY,
            self.variables["t_vec"],
            self.variables["deltaNM_noshots"],
        )

    def test_sem_uses_all_times_and_the_original_sample_count(self):
        sem = self.load_sem()
        self.assertEqual(sem.shape, (76,))
        np.testing.assert_allclose(
            sem, np.sqrt(self.summary["trajectory_variance_ddof1"] / 5_000_000),
            rtol=1e-12, atol=1e-14,
        )
        self.assertEqual(sem[0], 0)
        self.assertAlmostEqual(sem[-1], 0.0003597183257509643)

    def test_invalid_or_mismatched_data_are_rejected(self):
        cases = (
            ("time", 0.01, "time grid"),
            ("exact_mean", 0.0, "archived errors"),
            ("trajectory_variance_ddof1", -1.0, "Negative"),
            ("mean_se_extrapolated_N5000000", -1.0, "Negative"),
            ("mean_se_extrapolated_N5000000", 0.1, "sample count"),
            ("trajectory_variance_ddof1", np.nan, "Non-finite"),
        )
        for column, value, message in cases:
            with self.subTest(column=column, value=value):
                invalid = self.summary.copy()
                invalid[column][1] = value
                with patch.object(figures.np, "genfromtxt", return_value=invalid):
                    with self.assertRaisesRegex(ValueError, message):
                        self.load_sem()
        for invalid in (self.summary[:-1], self.summary[["time"]]):
            with patch.object(figures.np, "genfromtxt", return_value=invalid):
                with self.assertRaisesRegex(ValueError, "shape or columns"):
                    self.load_sem()

    def test_plot_preserves_points_and_maps_the_signed_interval(self):
        selected = np.r_[0:76:2, 75]
        sem = self.load_sem()[selected]
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(figures.plt, "close"):
                figures.draw_error_figure(self.variables, Path(directory) / "noise.png")
                figure = figures.plt.gcf()
            try:
                axes = figure.axes[0]
                ours = axes.containers[0]
                np.testing.assert_array_equal(
                    axes.lines[0].get_ydata(), self.variables["deltaJ0_noshots"][selected]
                )
                np.testing.assert_array_equal(
                    ours.lines[0].get_xdata(), self.variables["t_vec"][selected]
                )
                np.testing.assert_array_equal(
                    ours.lines[0].get_ydata(), self.variables["deltaNM_noshots"][selected]
                )
                deviation = (
                    self.summary["archive_mean_N5000000"] - self.summary["exact_mean"]
                )[selected]
                left, right = deviation - sem, deviation + sem
                lower = np.where(
                    left * right <= 0, 0, np.minimum(np.abs(left), np.abs(right))
                )
                upper = np.maximum(np.abs(left), np.abs(right))
                segments = np.asarray(ours.lines[2][0].get_segments())
                np.testing.assert_allclose(segments[:, 0, 1], lower, atol=1e-12)
                np.testing.assert_allclose(segments[:, 1, 1], upper, atol=1e-12)
                self.assertTrue(np.all(segments[:, 0, 1] >= 0))
                np.testing.assert_array_equal(figure.get_size_inches(), [8, 6])
                self.assertEqual(axes.get_xlim(), (0, 1.52))
                self.assertEqual(axes.get_ylim(), (0, 0.002))
                self.assertEqual(axes.get_ylabel(), "Observable estimation error")
                self.assertEqual(axes.yaxis.label.get_fontsize(), 20)
                self.assertEqual(ours.lines[0].get_markersize(), 10)
                self.assertEqual(ours.lines[0].get_color(), figures.ORANGE)
            finally:
                figures.plt.close(figure)


if __name__ == "__main__":
    unittest.main()
