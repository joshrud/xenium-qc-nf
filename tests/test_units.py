"""Fast unit tests (no Xenium data needed): transform fitting and threshold grading."""
import numpy as np

from xenium_qc.alignment import apply, describe, fit_similarity
from xenium_qc.run_metrics import grade


def test_fit_similarity_recovers_rotation_scale_shift():
    rng = np.random.default_rng(0)
    he = rng.uniform(0, 20000, (20, 2))
    th, s = np.radians(-90), 1.29
    M = np.array([[s * np.cos(th), -s * np.sin(th), -11900], [s * np.sin(th), s * np.cos(th), 28900], [0, 0, 1]])
    fit = fit_similarity(he, apply(M, he))
    assert np.allclose(fit, M, atol=1e-6)
    d = describe(fit)
    assert abs(d["scale"] - s) < 1e-9 and abs(d["rotation_deg"] + 90) < 1e-6 and not d["reflected"]


def test_grade_directions():
    assert grade(0.81, "min", 0.6, 0.5) == "pass"
    assert grade(0.55, "min", 0.6, 0.5) == "WARNING"
    assert grade(0.40, "min", 0.6, 0.5) == "ERROR"
    assert grade(0.03, "max", 0.025, 0.05) == "WARNING"
    assert grade(float("nan"), "max", 0.025, 0.05) == "n/a"
