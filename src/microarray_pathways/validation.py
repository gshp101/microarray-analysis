"""Check the drug rankings against known clinical trial outcomes.

The ranking routes never see trial data, so a benchmark of drugs with positive and
failed trials in one disease asks two things of each route:

- discrimination: do drugs with positive trials rank above drugs whose trials failed
  (AUC, Mann-Whitney), and above the rest of the drugs scored?
- specificity: does the disease's own dataset rank its positive drugs higher than the
  other datasets do? A route that ranks them high everywhere is picking up generic
  biology (or generic drug properties), not the disease signature.

Ranks are percentiles within each route and dataset, 100 = top of the list.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .treatments import _drug_key

BENCHMARKS = Path(__file__).parent / "benchmarks"

# route -> (file, name column, score column, lower score ranks higher)
ROUTES = {
    "kegg_targets": ("drugs_targets.tsv", "drug", None, None),
    "l1000_reversal": ("drugs_reversal.tsv", "compound", "z_metabolic", True),
    "proximity": ("drugs_proximity.tsv", "drug", "z", True),
}


def load_benchmark(path: str | Path) -> pd.DataFrame:
    """Benchmark table with 'drug' and 'outcome' (positive / negative) columns."""
    bench = pd.read_csv(path, sep="\t")
    bad = set(bench["outcome"]) - {"positive", "negative"}
    if bad:
        raise ValueError(f"outcome must be positive or negative, got {sorted(bad)}")
    bench["key"] = bench["drug"].map(_drug_key)
    return bench


def find_file(roots: list[Path], accession: str, name: str) -> Path | None:
    """First roots[i]/<accession>/<name> that exists (report output, or results_drugs/ etc.)."""
    for root in roots:
        path = Path(root) / accession / name
        if path.exists():
            return path
    return None


def route_percentiles(table: pd.DataFrame, route: str) -> pd.Series:
    """Percentile rank (100 = top) per normalized drug name; a name listed twice keeps its best row.

    KEGG target tables are already in ranked order (opposes change, marketed, score), so
    row position is the rank. The other routes rank by their score column.
    """
    _, name_col, score_col, ascending = ROUTES[route]
    if score_col is None:
        order = pd.Series(np.arange(len(table), 0, -1, dtype=float), index=table.index)
    else:
        score = table[score_col].astype(float)
        order = -score if ascending else score
    pct = order.rank(pct=True) * 100
    keys = table[name_col].map(_drug_key)
    return pct.groupby(keys).max()


def auc(high: pd.Series, low: pd.Series) -> tuple[float, float]:
    """P(a random 'high' drug outranks a random 'low' one) and its one-sided Mann-Whitney p."""
    if len(high) == 0 or len(low) == 0:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(high, low, alternative="greater")
    return u / (len(high) * len(low)), p


def validate(bench: pd.DataFrame, accessions: list[str], target: str,
             roots: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score the benchmark in every dataset; `target` is the dataset of the benchmark's disease.

    Returns (ranks, summary). ranks: one row per benchmark drug, route and dataset.
    summary: per route, coverage, AUCs in the target dataset, and how far the target
    dataset lifts positive drugs over the median of the other datasets.
    """
    rows, tables = [], {}
    for acc in accessions:
        for route, (fname, *_) in ROUTES.items():
            path = find_file(roots, acc, fname)
            if path is None:
                continue
            pct = route_percentiles(pd.read_csv(path, sep="\t"), route)
            tables[acc, route] = pct
            for _, b in bench.iterrows():
                rows.append({"dataset": acc, "route": route, "drug": b["drug"], "outcome": b["outcome"],
                             "scored": b["key"] in pct.index, "percentile": pct.get(b["key"], np.nan),
                             "n_scored": len(pct)})
    ranks = pd.DataFrame(rows)
    if ranks.empty:
        raise FileNotFoundError(f"no drug tables found for {accessions} under {roots}")

    summary = []
    for route, grp in ranks.groupby("route", sort=False):
        here = grp[(grp["dataset"] == target) & grp["scored"]]
        if here.empty:
            continue
        pos = here.loc[here["outcome"] == "positive", "percentile"]
        neg = here.loc[here["outcome"] == "negative", "percentile"]
        auc_pn, p_pn = auc(pos, neg)
        auc_pr, p_pr = auc(pos, tables[target, route].drop(bench["key"], errors="ignore"))
        other = grp[(grp["dataset"] != target) & grp["scored"]]
        lift = (here.set_index("drug")["percentile"]
                - other.groupby("drug")["percentile"].median()).dropna()
        lifts = {o: lift.reindex(here.loc[here["outcome"] == o, "drug"]).dropna()
                 for o in ("positive", "negative")}
        summary.append({
            "route": route, "n_scored": int(here["n_scored"].iloc[0]),
            "positives_scored": f"{len(pos)}/{(bench['outcome'] == 'positive').sum()}",
            "negatives_scored": f"{len(neg)}/{(bench['outcome'] == 'negative').sum()}",
            "median_pct_positive": pos.median(), "median_pct_negative": neg.median(),
            "auc_pos_vs_neg": auc_pn, "p_pos_vs_neg": p_pn,
            "auc_pos_vs_unlabeled": auc_pr, "p_pos_vs_unlabeled": p_pr,
            "other_datasets": other["dataset"].nunique(),
        })
        for o, lo in lifts.items():
            summary[-1][f"median_lift_{o}"] = lo.median() if len(lo) else np.nan
            summary[-1][f"p_lift_{o}"] = (stats.wilcoxon(lo, alternative="greater").pvalue
                                          if len(lo) >= 5 else np.nan)
    return ranks, pd.DataFrame(summary)


def write(ranks: pd.DataFrame, summary: pd.DataFrame, target: str, out_dir: str | Path) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    wide = ranks.pivot_table(index=["drug", "outcome", "route"], columns="dataset",
                             values="percentile", dropna=False).round(1)
    wide = wide[[target] + [c for c in wide.columns if c != target]].reset_index()
    wide.sort_values(["route", "outcome", target], ascending=[True, False, False]).to_csv(
        out / "benchmark_ranks.tsv", sep="\t", index=False)
    summary.round(4).to_csv(out / "benchmark_summary.tsv", sep="\t", index=False)
    return out
