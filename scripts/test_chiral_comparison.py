"""Checks for shared circular scoring and exact Figure 1 input matching."""

import unittest

import numpy as np
import pandas as pd

from supplementary_chiral_comparison import matched_panel, score_phase


class ComparableScoringTests(unittest.TestCase):
    def test_rotation_reflection_and_wrap_invariance(self):
        ct = np.arange(0, 24, 3, dtype=float)
        truth = ct * 2 * np.pi / 24
        for phase in (truth + 5.7, -truth + 5.7, truth + 8 * np.pi):
            metrics, aligned, error = score_phase(phase, ct)
            self.assertAlmostEqual(metrics["R"], 1)
            self.assertAlmostEqual(metrics["mean_abs_err_h"], 0)
            np.testing.assert_allclose(error, 0, atol=1e-12)
            np.testing.assert_allclose(np.exp(1j * aligned), np.exp(1j * truth), atol=1e-12)

    def test_known_symmetric_phase_error(self):
        ct = np.arange(0, 24, 3, dtype=float)
        phase = ct * 2 * np.pi / 24 + np.tile([0.2, -0.2], 4)
        metrics, _, error = score_phase(phase, ct)
        self.assertAlmostEqual(metrics["R"], np.cos(0.2))
        self.assertAlmostEqual(metrics["mean_abs_err_h"], 0.2 * 24 / (2 * np.pi))
        np.testing.assert_allclose(np.abs(error), metrics["mean_abs_err_h"])

    def test_invalid_phases_rejected(self):
        for phase in ([0, 1], [0, np.nan, 2], [0, np.inf, 2]):
            with self.assertRaises(ValueError):
                score_phase(phase, [0, 8, 16])

    def test_panel_reconstruction_and_mismatch_detection(self):
        ct = np.arange(0, 24, 3)
        names = [f"Adr_CT{value}" for value in ct]
        ids = [str(value) for value in range(6)]
        values = np.exp(np.random.default_rng(17).normal(size=(6, 8)))
        expr = pd.DataFrame(values, index=ids, columns=names)
        angles = 2 * np.pi * np.arange(6) / 6
        ft = pd.DataFrame({"ID_REF": ids, "acs_position": np.arange(6), "acs_angle": angles})
        clr = np.log(values).T
        clr -= clr.mean(axis=1, keepdims=True)
        clr -= clr.mean(axis=0, keepdims=True)
        xy = clr @ np.array([np.cos(angles), np.sin(angles)]).T * np.sqrt(2 / 6)
        st = pd.DataFrame({"sample": names, "ct": ct, "acs_u": xy[:, 0], "acs_v": xy[:, 1]})
        matrix, error = matched_panel(expr, ft.iloc[::-1], st)
        self.assertLess(error, 1e-12)
        np.testing.assert_allclose(matrix, np.log2(values))
        altered = expr.copy()
        altered.iloc[0, 0] *= 2
        with self.assertRaisesRegex(ValueError, "reproduce saved"):
            matched_panel(altered, ft, st)
        st.loc[0, "ct"] = 2
        with self.assertRaisesRegex(ValueError, "collection times"):
            matched_panel(expr, ft, st)


if __name__ == "__main__":
    unittest.main()
