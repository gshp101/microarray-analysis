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


ROUTE_LABELS = {"kegg_targets": "KEGG drug targets", "l1000_reversal": "L1000 signature reversal",
                "proximity": "Network proximity", "any_route": "Best of the three routes"}

# trial-outcome benchmark for each curated dataset's disease
DEFAULT_BENCHMARKS = {
    "GSE14520": "hcc_trials.tsv",        # liver cancer (HCC)
    "GSE10072": "nsclc_trials.tsv",      # lung adenocarcinoma
    "GSE15852": "breast_trials.tsv",     # breast cancer
    "GSE32323": "crc_trials.tsv",        # colorectal cancer
    "GSE5281": "alzheimers_trials.tsv",  # Alzheimer's disease
    "GSE48350": "alzheimers_trials.tsv",
}


def benchmark_scores(bench: pd.DataFrame, accession: str, roots: list[Path]) -> pd.DataFrame:
    """Percentile of every benchmark drug in each route of one dataset, plus the best of the routes.

    Every benchmark drug is kept: one a route did not score was not identified by it, so
    it gets percentile 0 (tied at the bottom).
    """
    rows = []
    for route, (fname, *_) in ROUTES.items():
        path = find_file(roots, accession, fname)
        if path is None:
            continue
        pct = route_percentiles(pd.read_csv(path, sep="\t"), route)
        for _, b in bench.iterrows():
            rows.append({"dataset": accession, "route": route, "drug": b["drug"], "outcome": b["outcome"],
                         "scored": b["key"] in pct.index, "percentile": float(pct.get(b["key"], 0.0))})
    if not rows:
        return pd.DataFrame(rows)
    out = pd.DataFrame(rows)
    best = (out.groupby(["dataset", "drug", "outcome"], as_index=False, sort=False)
            .agg(scored=("scored", "any"), percentile=("percentile", "max")))
    best["route"] = "any_route"
    return pd.concat([out, best[out.columns]], ignore_index=True)


def roc_table(sets: list[tuple[str, pd.DataFrame]], roots: list[Path],
              threshold: float = 90.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ROC of each route's ranking against trial outcomes, pooled over (accession, benchmark) sets.

    Each dataset is scored against its own disease's benchmark; the (drug, dataset) pairs
    are then pooled. Detections and false alarms count positive and negative pairs at or
    above `threshold` (percentile, 100 = top).

    Returns (scores, summary): scores has one row per drug, dataset and route; summary per
    route has AUC, detections, false alarms and the ROC points (fpr, tpr lists).
    """
    parts = [benchmark_scores(b, acc, roots) for acc, b in sets]
    parts = [p for p in parts if not p.empty]
    if not parts:
        raise FileNotFoundError(f"no drug tables found for {[a for a, _ in sets]} under {roots}")
    scores = pd.concat(parts, ignore_index=True)
    scores["flagged"] = scores["percentile"] >= threshold
    summary = []
    for route, grp in scores.groupby("route", sort=False):
        score = grp["percentile"].to_numpy(dtype=float)
        truth = (grp["outcome"] == "positive").to_numpy()
        scored, flagged = grp["scored"].to_numpy(), grp["flagged"].to_numpy()
        fpr, tpr = roc_points(score, truth)
        a, p = auc(pd.Series(score[truth]), pd.Series(score[~truth]))
        summary.append({"route": route, "datasets": grp["dataset"].nunique(),
                        "n_positive": int(truth.sum()), "n_negative": int((~truth).sum()),
                        "positives_scored": int((scored & truth).sum()),
                        "negatives_scored": int((scored & ~truth).sum()),
                        "auc": a, "p": p, "threshold_pct": threshold,
                        "detections": int((flagged & truth).sum()),
                        "false_alarms": int((flagged & ~truth).sum()),
                        "fpr": fpr, "tpr": tpr})
    return scores, pd.DataFrame(summary)


def roc_points(score: np.ndarray, truth: np.ndarray) -> tuple[list[float], list[float]]:
    """ROC curve through every distinct score, ties stepping diagonally (start at 0,0)."""
    pos, neg = truth.sum(), (~truth).sum()
    fpr, tpr = [0.0], [0.0]
    for cut in np.unique(score)[::-1]:
        hit = score >= cut
        fpr.append((hit & ~truth).sum() / neg)
        tpr.append((hit & truth).sum() / pos)
    return fpr, tpr


def _n(k: int, word: str) -> str:
    return f"{k} {word}" + ("" if k == 1 else "s")


def plot_roc(summary: pd.DataFrame, title: str, path: str | Path, subtitle: str = "") -> Path:
    from . import plots  # matplotlib only when plotting

    unit = "drugs" if summary["datasets"].max() == 1 else "drug-dataset pairs"

    plt = plots.plt
    n = len(summary)
    cols = 2 if n == 4 else n
    rows = -(-n // cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 4.9 * rows + 0.4), dpi=150, squeeze=False,
                             layout="constrained")
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    for ax, (_, r) in zip(axes.flat, summary.iterrows()):
        plots._style(ax)
        ax.grid(axis="both", color="#e6e5e0", linewidth=0.8)
        ax.plot([0, 1], [0, 1], color=plots.NEUTRAL, linewidth=1, linestyle="--")
        ax.plot(r["fpr"], r["tpr"], color=plots.DOWN, linewidth=2)
        # operating point at the threshold
        x, y = r["false_alarms"] / r["n_negative"], r["detections"] / r["n_positive"]
        ax.plot([x], [y], "o", color=plots.UP, markersize=6, zorder=3)
        ax.annotate(f"top {100 - r['threshold_pct']:.0f}%", (x, y), xytext=(6, -12),
                    textcoords="offset points", fontsize=8, color=plots.MUTED)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect("equal")
        ax.set_xlabel(f"False alarm rate (of {r['n_negative']} failed-trial {unit})", fontsize=9, color=plots.MUTED)
        ax.set_ylabel(f"Detection rate (of {r['n_positive']} positive-trial {unit})", fontsize=9, color=plots.MUTED)
        ax.set_title(f"{ROUTE_LABELS.get(r['route'], r['route'])}\n"
                     f"AUC {r['auc']:.2f} · {_n(r['detections'], 'detection')} · {_n(r['false_alarms'], 'false alarm')}",
                     loc="left", fontsize=10, color=plots.INK)
    fig.suptitle(title + (f"\n{subtitle}" if subtitle else ""), x=0.01, ha="left", fontsize=11, color=plots.INK)
    fig.savefig(path)
    plt.close(fig)
    return Path(path)


def write_roc(scores: pd.DataFrame, summary: pd.DataFrame, title: str, out_dir: str | Path) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    thr = summary["threshold_pct"].iloc[0]
    plot_roc(summary, title, out / "roc.png",
             subtitle=f"Counts at the top {100 - thr:.0f}% of each ranking; drugs a route did not score rank last.")
    scores.sort_values(["route", "dataset", "percentile"], ascending=[True, True, False]).to_csv(
        out / "roc_scores.tsv", sep="\t", index=False)
    summary.drop(columns=["fpr", "tpr"]).round(4).to_csv(out / "roc_summary.tsv", sep="\t", index=False)
    return out
