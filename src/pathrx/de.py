"""Differential expression between two sample groups.

For each gene: log2 fold change (case mean - control mean on the log2 scale),
Welch's t-test, and Benjamini-Hochberg false discovery rate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .geo import ExpressionDataset


def differential_expression(
    expression: pd.DataFrame,
    groups: pd.Series,
    case: str,
    control: str,
    min_mean: float | None = None,
) -> pd.DataFrame:
    """Compare `case` vs `control` samples gene by gene.

    expression: genes x samples, log2 scale.
    groups: sample -> group label (samples with other labels or NaN are ignored).
    min_mean: optionally drop genes whose mean log2 expression across both
        groups is below this value (low-expressed genes are mostly noise).

    Returns one row per gene sorted by p-value, with columns
    log2fc, mean_case, mean_control, t, pvalue, fdr.
    """
    groups = groups.reindex(expression.columns)
    a = expression.loc[:, (groups == case).to_numpy()]
    b = expression.loc[:, (groups == control).to_numpy()]
    if a.shape[1] < 2 or b.shape[1] < 2:
        raise ValueError(f"need at least 2 samples per group, got {case}={a.shape[1]}, "
                         f"{control}={b.shape[1]}")

    mean_a, mean_b = a.mean(axis=1), b.mean(axis=1)
    keep = pd.Series(True, index=expression.index)
    if min_mean is not None:
        keep &= (mean_a + mean_b) / 2 >= min_mean
    keep &= (a.notna().sum(axis=1) >= 2) & (b.notna().sum(axis=1) >= 2)
    a, b, mean_a, mean_b = a[keep], b[keep], mean_a[keep], mean_b[keep]

    t, p = stats.ttest_ind(a.to_numpy(), b.to_numpy(), axis=1, equal_var=False, nan_policy="omit")
    t, p = np.asarray(t, dtype=float), np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    fdr = np.full_like(p, np.nan)
    fdr[ok] = multipletests(p[ok], method="fdr_bh")[1]

    res = pd.DataFrame({
        "log2fc": mean_a - mean_b,
        "mean_case": mean_a,
        "mean_control": mean_b,
        "t": t,
        "pvalue": p,
        "fdr": fdr,
    }, index=a.index)
    res.attrs.update(case=case, control=control, n_case=a.shape[1], n_control=b.shape[1])
    return res.sort_values("pvalue")


def run(ds: ExpressionDataset, case: str, control: str, min_mean: float | None = None) -> pd.DataFrame:
    return differential_expression(ds.expression, ds.samples["group"], case, control, min_mean)


def significant(de: pd.DataFrame, fdr: float = 0.05, min_abs_lfc: float = 1.0) -> pd.DataFrame:
    return de[(de["fdr"] < fdr) & (de["log2fc"].abs() >= min_abs_lfc)]
