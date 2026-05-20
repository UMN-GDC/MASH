"""AdjHE estimator — adjusted Haseman-Elston heritability estimation."""

import logging
import numpy as np
from scipy.linalg import block_diag
import pandas as pd


def AdjHE(A, df, mp, random_groups=None, npc=0, std=False):
    A = np.asarray(A, dtype=float)

    if random_groups is None:
        sigma_g, sigma_e, n, trA, trA2 = _adjhe_2comp(A, df, mp, npc, std)
    else:
        df = df.sort_values(random_groups).dropna(subset=[random_groups])
        A = A[df.index][:, df.index]
        sigma_g, sigma_e, n, trA, trA2 = _adjhe_3comp(A, df, mp, random_groups, std)

    var_h2 = _h2_variance(n, trA, trA2)
    return _package(sigma_g, sigma_e, var_h2)


def _extract_y(df, mp, std):
    y = np.array(df[mp], dtype=float)
    if std:
        y = (y - y.mean()) / y.std()
    return y


def _q_projections(df, npc, A, y, n):
    """Build Q = I - V(V'V)^{-1}V' via QR and compute projected traces without materializing Q."""
    pc_cols = [c for c in df.columns if c.startswith("pc")]
    V = np.array(df[pc_cols])[:, :npc]
    V, _ = np.linalg.qr(V, mode="reduced")

    Ay = A @ y
    AV = A @ V

    yV = y @ V
    VAy = V.T @ Ay
    N = V.T @ AV
    N2 = AV.T @ AV

    trQA = np.trace(A) - np.trace(N)
    trQAQA = np.trace(A @ A) - 2 * np.trace(N2) + np.trace(N @ N)
    yQy = y @ y - yV @ yV
    yQAQy = y @ Ay - 2 * yV @ VAy + yV @ N @ yV

    return trQA, trQAQA, yQy, yQAQy, n - npc


def _adjhe_2comp(A, df, mp, npc, std):
    y = _extract_y(df, mp, std)
    n = A.shape[0]
    trA = np.trace(A)
    trA2 = np.trace(A @ A)
    yAy = y @ A @ y
    yty = y @ y

    if npc > 0:
        trA_adj, trA2_adj, yty_adj, yAy_adj, n_adj = _q_projections(df, npc, A, y, n)
    else:
        trA2_adj, trA_adj, n_adj, yAy_adj, yty_adj = trA2, trA, n, yAy, yty

    det = trA2_adj * n_adj - trA_adj ** 2

    if det <= 0 or np.isnan(det) or np.isinf(det):
        logging.warning(
            f"AdjHE: Non-positive or invalid determinant ({det:.6e}) with npc={npc}. "
            f"trA2={trA2:.4e}, trA={trA:.4e}, n={n}, "
            f"trA2_adj={trA2_adj:.4e}, trA_adj={trA_adj:.4e}, n_adj={n_adj}, "
            f"yAy={yAy:.4e}, yty={yty:.4e}, "
            f"yAy_adj={yAy_adj:.4e}, yty_adj={yty_adj:.4e}"
        )
        return np.nan, np.nan, n, trA, trA2

    sigma_g = (yAy_adj * n_adj - trA_adj * yty_adj) / det
    sigma_e = (trA2_adj * yty_adj - trA_adj * yAy_adj) / det

    if sigma_g < 0:
        logging.warning(
            f"AdjHE: Negative genetic variance. "
            f"sigmas=({sigma_g:.4e}, {sigma_e:.4e}), "
            f"yAy={yAy:.4e}, yty={yty:.4e}, trA2={trA2:.4e}, trA={trA:.4e}, "
            f"topleft={trA2_adj:.4e}, offdiag={trA_adj:.4e}, bottomright={n_adj:.4e}"
        )

    return sigma_g, sigma_e, n, trA, trA2


def _adjhe_3comp(A, df, mp, random_groups, std):
    y = _extract_y(df, mp, std)
    n = A.shape[0]

    proj_cols = [c for c in df.columns if c.startswith("pc")]
    proj_cols.append(random_groups)
    X = np.array(pd.get_dummies(df[proj_cols]))

    q, _ = np.linalg.qr(X)
    Q = np.eye(n) - q @ q.T

    sizes = np.unique(df[random_groups], return_counts=True)[1]
    S = block_diag(*[np.ones((s, s)) for s in sizes])

    QSQ = Q @ S @ Q
    youter = np.outer(y, y)

    trA = np.trace(A)
    trA2 = np.trace(A @ A)
    trQSQA = np.trace(QSQ @ A)
    trQSQ = np.trace(QSQ)
    trQSQQSQ = np.trace(QSQ @ QSQ)

    XtX = np.array([
        [trA2, trQSQA, trA],
        [trQSQA, trQSQQSQ, trQSQ],
        [trA, trQSQ, n],
    ])
    Xty = np.array([
        np.trace(A @ youter),
        np.trace(QSQ @ youter),
        np.trace(youter),
    ])

    try:
        sigmas = np.linalg.solve(XtX, Xty)
    except np.linalg.LinAlgError:
        logging.warning("AdjHE (random_groups): Singular system in 3-component solve")
        return np.nan, np.nan, n, trA, trA2

    if np.any(np.isnan(sigmas)):
        logging.warning("AdjHE (random_groups): NaN in solve result")
        return np.nan, np.nan, n, trA, trA2

    if sigmas[0] < 0:
        logging.warning(
            f"AdjHE (random_groups): Negative genetic variance. "
            f"sigmas=({sigmas[0]:.4e}, {sigmas[1]:.4e}, {sigmas[2]:.4e})"
        )

    return sigmas[0], sigmas[2], n, trA, trA2


def _h2_variance(n, trA, trA2):
    v = 2 * n / (n * trA2 - trA ** 2)
    if v <= 0:
        logging.warning(
            f"AdjHE: Variance estimate is negative ({v:.4e}), setting as absolute value"
        )
        v = abs(v)
    return v


def _package(sigma_g, sigma_e, var_h2):
    if np.isnan(sigma_g):
        logging.warning("AdjHE: Estimate is NaN (singular system, returning NaN)")
        return {"h2": np.nan, "var(h2)": var_h2}

    h2 = sigma_g / (sigma_g + sigma_e)

    if h2 < 0:
        logging.warning(
            f"AdjHE: h2 clamped from {h2:.4e} to 0. "
            "This indicates heritability near zero or a poor model fit"
        )
        h2 = 0
    elif h2 > 1:
        logging.warning(f"AdjHE: h2 clamped from {h2:.4e} to 1")
        h2 = 1

    return {"h2": h2, "var(h2)": var_h2}
