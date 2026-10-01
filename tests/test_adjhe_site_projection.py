"""Regression test: site indicators must not enter the PC projection.

In the 3-component AdjHE path, S is built from the site indicator matrix D
(S = block_diag of ones == D D'). The original code also appended random_groups
to proj_cols, so q spanned the columns of D, Q D = 0, and

    QSQ = Q D D' Q = (Q D)(Q D)' = 0

identically. XtX then lost its middle row, cond() = inf, and every single RE
phenotype came back flagged ill_conditioned at every npc and sample size
(2080/2080 for gordon, 400/400 for probaConns, 36/36 and 36/36 for SA).

A component being estimated as a random effect must not also be projected away
as a fixed effect. This test pins both halves of the fix: the projection no
longer contains the site dummies, and npc actually selects columns.
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.linalg import block_diag

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

from Estimate.estimators.AdjHE import AdjHE  # noqa: E402


def _make_data(n=120, n_sites=8, n_pc=20, seed=0):
    rng = np.random.default_rng(seed)
    M = rng.normal(size=(n, n))
    A = (M + M.T) / 2.0
    A = A / np.sqrt(np.einsum("ij,ij->i", A, A))[:, None]
    np.fill_diagonal(A, 1.0)

    site = np.repeat(np.arange(n_sites), n // n_sites)
    df = pd.DataFrame({f"pc{i + 1}": rng.normal(size=n) for i in range(n_pc)})
    df["site"] = [f"s{v}" for v in site]
    df["pheno"] = A @ rng.normal(size=n) + 0.5 * site + rng.normal(size=n)
    return A, df


def test_site_in_projection_annihilates_qsq():
    """The original construction is exactly degenerate - documents the bug."""
    _, df = _make_data()
    n = len(df)

    proj_cols = [c for c in df.columns if c.startswith("pc")]
    proj_cols.append("site")
    X = np.array(pd.get_dummies(df[proj_cols]))
    q, _ = np.linalg.qr(X)
    Q = np.eye(n) - q @ q.T

    sizes = np.unique(df["site"], return_counts=True)[1]
    S = block_diag(*[np.ones((s, s)) for s in sizes])
    QSQ = Q @ S @ Q

    assert np.abs(QSQ).max() < 1e-8, "site in proj_cols must make QSQ vanish"


def test_random_groups_yields_finite_estimate():
    """RE no longer returns all-NaN with flag ill_conditioned."""
    A, df = _make_data()
    out = AdjHE(A=A, df=df, mp="pheno", random_groups="site", npc=10, std=False)

    assert out["flag"] != "ill_conditioned"
    assert np.isfinite(out["G"]), "genetic variance must be finite"
    assert np.isfinite(out["S"]), "site variance must be estimable"
    assert np.isfinite(out["E"]), "residual variance must be finite"


def test_npc_selects_columns():
    """npc must change the projection, otherwise npc=10 and npc=20 are copies."""
    A, df = _make_data()
    low = AdjHE(A=A, df=df, mp="pheno", random_groups="site", npc=5, std=False)
    high = AdjHE(A=A, df=df, mp="pheno", random_groups="site", npc=20, std=False)

    assert low["flag"] != "ill_conditioned" and high["flag"] != "ill_conditioned"
    assert not np.isclose(low["G"], high["G"]), "npc=5 and npc=20 must differ"


def test_npc_zero_does_not_crash():
    A, df = _make_data()
    out = AdjHE(A=A, df=df, mp="pheno", random_groups="site", npc=0, std=False)
    assert out["flag"] != "ill_conditioned"
    assert np.isfinite(out["G"])


def test_two_component_path_untouched():
    """FE (random_groups=None) must behave as before."""
    A, df = _make_data()
    out = AdjHE(A=A, df=df, mp="pheno", random_groups=None, npc=10, std=False)
    assert out["flag"] != "ill_conditioned"
    assert np.isfinite(out["h2"])
    assert 0.0 <= out["h2"] <= 1.0
