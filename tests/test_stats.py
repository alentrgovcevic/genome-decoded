"""The small statistical pieces the report's numbers rest on."""

import numpy as np

from ancestry_membership import lda, posterior
from pgs import pct


def test_pct():
    dist = np.arange(100.0)
    assert pct(-1, dist) == 0
    assert pct(50, dist) == 50
    assert pct(1000, dist) == 100


def test_lda_assigns_to_nearest_cluster():
    rng = np.random.default_rng(1)
    X = np.vstack([rng.normal(0, 1, (50, 2)), rng.normal(6, 1, (50, 2))])
    y = np.array(["a"] * 50 + ["b"] * 50)
    means, icov = lda(X, y, ["a", "b"])
    p, m2 = posterior(np.array([5.8, 6.1]), means, icov)
    assert p[1] > 0.99 and np.isclose(p.sum(), 1)
    p, _ = posterior(means.mean(axis=0), means, icov)  # equidistant, equal priors: exactly 50/50
    assert np.isclose(p[0], 0.5)
