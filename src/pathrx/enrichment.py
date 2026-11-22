"""Score pathways for a differential expression result.

Two standard methods:
  ora   over-representation: are the significantly up (or down) genes
        enriched in a pathway, versus all measured genes? (hypergeometric test)
  gsea  ranked gene set enrichment on every measured gene, ranked by t
        statistic (gseapy prerank). Uses all genes, so it picks up pathways
        where many genes shift a little, which is typical for metabolism.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _restrict(gene_sets: dict[str, set[str]], universe: set[str], min_size: int, max_size: int):
    out = {}
    for name, genes in gene_sets.items():
        g = genes & universe
        if min_size <= len(g) <= max_size:
            out[name] = g
    return out


def ora(
    selected: set[str],
    universe: set[str],
    gene_sets: dict[str, set[str]],
    min_size: int = 10,
    max_size: int = 500,
) -> pd.DataFrame:
    """Hypergeometric over-representation of `selected` genes in each gene set.

    The background is genes that were measured AND belong to at least one gene
    set, the usual choice for pathway ORA.
    """
    sets = _restrict(gene_sets, universe, min_size, max_size)
    background = set().union(*sets.values()) if sets else set()
    hits_all = selected & background
    N, n = len(background), len(hits_all)
    rows = []
    for name, genes in sets.items():
        hits = genes & hits_all
        k, K = len(hits), len(genes)
        p = stats.hypergeom.sf(k - 1, N, K, n) if k else 1.0
        rows.append({"pathway": name, "size": K, "hits": k,
                     "expected": K * n / N if N else np.nan,
                     "pvalue": p, "genes": ",".join(sorted(hits))})
    res = pd.DataFrame(rows, columns=["pathway", "size", "hits", "expected", "pvalue", "genes"])
    if len(res):
        res["fdr"] = multipletests(res["pvalue"], method="fdr_bh")[1]
    return res.sort_values("pvalue").reset_index(drop=True)


def ora_up_down(de: pd.DataFrame, gene_sets: dict[str, set[str]], fdr: float = 0.05,
                min_abs_lfc: float = 1.0, **kw) -> pd.DataFrame:
    """ORA run separately on up- and down-regulated genes; one table with a direction column."""
    sig = (de["fdr"] < fdr) & (de["log2fc"].abs() >= min_abs_lfc)
    universe = set(de.index)
    parts = []
    for direction, mask in (("up", sig & (de["log2fc"] > 0)), ("down", sig & (de["log2fc"] < 0))):
        r = ora(set(de.index[mask]), universe, gene_sets, **kw)
        r.insert(1, "direction", direction)
        parts.append(r)
    return pd.concat(parts).sort_values("pvalue").reset_index(drop=True)


def gsea(
    de: pd.DataFrame,
    gene_sets: dict[str, set[str]],
    rank_by: str = "t",
    min_size: int = 10,
    max_size: int = 500,
    permutations: int = 1000,
    seed: int = 0,
    threads: int = 4,
) -> pd.DataFrame:
    """Preranked GSEA. Positive NES means the pathway is shifted up in the case group."""
    import gseapy

    rnk = de[rank_by].dropna().sort_values(ascending=False)
    sets = {k: sorted(v) for k, v in _restrict(gene_sets, set(rnk.index), min_size, max_size).items()}
    pre = gseapy.prerank(rnk=rnk, gene_sets=sets, min_size=min_size, max_size=max_size,
                         permutation_num=permutations, seed=seed, threads=threads,
                         outdir=None, verbose=False)
    r = pre.res2d.rename(columns={"Term": "pathway", "NES": "nes", "ES": "es",
                                  "NOM p-val": "pvalue", "FDR q-val": "fdr",
                                  "Lead_genes": "leading_edge", "Tag %": "tag_pct"})
    for c in ("nes", "es", "pvalue", "fdr"):
        r[c] = pd.to_numeric(r[c], errors="coerce")
    r["size"] = r["pathway"].map(lambda p: len(sets[p]))
    r["direction"] = np.where(r["nes"] > 0, "up", "down")
    cols = ["pathway", "direction", "nes", "es", "pvalue", "fdr", "size", "leading_edge"]
    return r[cols].sort_values(["fdr", "pvalue"]).reset_index(drop=True)
