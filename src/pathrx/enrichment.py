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
    engine: str = "auto",
) -> pd.DataFrame:
    """Preranked GSEA. Positive NES means the pathway is shifted up in the case group.

    Uses gseapy when it is installed; otherwise (engine="numpy", or in the browser
    under Pyodide, where gseapy has no build) the NumPy implementation prerank_numpy,
    which follows the same method and gives close but not identical numbers.
    """
    rnk = de[rank_by].dropna().sort_values(ascending=False)
    sets = {k: sorted(v) for k, v in _restrict(gene_sets, set(rnk.index), min_size, max_size).items()}
    if engine == "auto":
        try:
            import gseapy  # noqa: F401
            engine = "gseapy"
        except ImportError:
            engine = "numpy"
    if engine == "numpy":
        r = prerank_numpy(rnk, sets, permutations=permutations, seed=seed)
    else:
        r = _prerank_gseapy(rnk, sets, min_size, max_size, permutations, seed, threads)
    r["size"] = r["pathway"].map(lambda p: len(sets[p]))
    r["direction"] = np.where(r["nes"] > 0, "up", "down")
    cols = ["pathway", "direction", "nes", "es", "pvalue", "fdr", "size", "leading_edge"]
    return r[cols].sort_values(["fdr", "pvalue"]).reset_index(drop=True)


def _prerank_gseapy(rnk, sets, min_size, max_size, permutations, seed, threads) -> pd.DataFrame:
    import gseapy

    pre = gseapy.prerank(rnk=rnk, gene_sets=sets, min_size=min_size, max_size=max_size,
                         permutation_num=permutations, seed=seed, threads=threads,
                         outdir=None, verbose=False)
    r = pre.res2d.rename(columns={"Term": "pathway", "NES": "nes", "ES": "es",
                                  "NOM p-val": "pvalue", "FDR q-val": "fdr",
                                  "Lead_genes": "leading_edge", "Tag %": "tag_pct"})
    for c in ("nes", "es", "pvalue", "fdr"):
        r[c] = pd.to_numeric(r[c], errors="coerce")
    return r


def _enrichment_scores(pos: np.ndarray, w: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Weighted Kolmogorov-Smirnov enrichment score for hit positions (rows of `pos`, sorted).

    pos: (m, k) positions of a set's genes in the ranked list; w: (m, k) their |rank metric|.
    The running sum peaks just after a hit and dips just before one, so only those
    points are evaluated. Returns (es, index of the hit where the peak or dip sits).
    """
    k = pos.shape[1]
    cw = np.cumsum(w, axis=1)
    total = cw[:, -1:]
    total = np.where(total > 0, total, 1.0)
    misses = (pos - np.arange(k)) / (n - k)
    after = cw / total - misses
    before = (cw - w) / total - misses
    i_max, i_min = after.argmax(axis=1), before.argmin(axis=1)
    rows = np.arange(len(pos))
    hi, lo = after[rows, i_max], before[rows, i_min]
    up = hi >= -lo
    return np.where(up, hi, lo), np.where(up, i_max, i_min)


def prerank_numpy(rnk: pd.Series, sets: dict[str, list[str]], permutations: int = 1000,
                  seed: int = 0, weight: float = 1.0, chunk: int = 100) -> pd.DataFrame:
    """Preranked GSEA (Subramanian et al. 2005) in NumPy, as gseapy.prerank computes it.

    rnk: gene -> rank metric, sorted descending. Null: gene labels permuted
    `permutations` times (each permutation shared by all sets, as in gseapy). NES divides
    ES by the mean null ES of the same sign; the nominal p-value is the share of same-sign
    null ES at least as extreme; FDR compares each NES with the pooled null NES of all
    sets (Subramanian's FDR q-value). Leading edge: the set's genes up to the peak (positive
    ES) or from the dip on (negative ES), as gseapy reports them.
    """
    genes = rnk.index.to_numpy()
    n = len(genes)
    weights = np.abs(rnk.to_numpy(dtype=float)) ** weight
    index = {g: i for i, g in enumerate(genes)}
    names = list(sets)
    members = [np.array(sorted(index[g] for g in sets[p] if g in index), dtype=np.int64) for p in names]

    rng = np.random.default_rng(seed)
    es = np.empty(len(names))
    leading = []
    null = np.empty((len(names), permutations))
    for j, idx in enumerate(members):
        e, at = _enrichment_scores(idx[None, :], weights[idx][None, :], n)
        es[j] = e[0]
        hit_genes = genes[idx]
        leading.append(";".join(hit_genes[:at[0] + 1] if e[0] >= 0 else hit_genes[at[0]:]))
    done = 0
    while done < permutations:
        m = min(chunk, permutations - done)
        perm = np.argsort(rng.random((m, n)), axis=1)  # perm[p, i]: rank of gene i in permutation p
        for j, idx in enumerate(members):
            pos = np.sort(perm[:, idx], axis=1)
            null[j, done:done + m] = _enrichment_scores(pos, weights[pos], n)[0]
        done += m

    pos_null, neg_null = np.where(null >= 0, null, np.nan), np.where(null < 0, null, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_pos = np.nanmean(pos_null, axis=1)
        mean_neg = -np.nanmean(neg_null, axis=1)
        nes = np.where(es >= 0, es / mean_pos, es / mean_neg)
        null_nes = np.where(null >= 0, null / mean_pos[:, None], null / mean_neg[:, None])
        pval = np.where(es >= 0,
                        (null >= es[:, None]).sum(axis=1) / (null >= 0).sum(axis=1),
                        (null <= es[:, None]).sum(axis=1) / (null < 0).sum(axis=1))

    flat = null_nes[np.isfinite(null_nes)]
    null_pos, null_neg = np.sort(flat[flat >= 0]), np.sort(flat[flat < 0])
    obs_pos, obs_neg = np.sort(nes[nes >= 0]), np.sort(nes[nes < 0])
    fdr = np.full(len(nes), np.nan)
    for j, v in enumerate(nes):
        if not np.isfinite(v):
            continue
        if v >= 0:
            a = (len(null_pos) - np.searchsorted(null_pos, v, "left")) / max(len(null_pos), 1)
            b = (len(obs_pos) - np.searchsorted(obs_pos, v, "left")) / max(len(obs_pos), 1)
        else:
            a = np.searchsorted(null_neg, v, "right") / max(len(null_neg), 1)
            b = np.searchsorted(obs_neg, v, "right") / max(len(obs_neg), 1)
        fdr[j] = min(a / b, 1.0) if b > 0 else 1.0
    return pd.DataFrame({"pathway": names, "nes": nes, "es": es, "pvalue": pval, "fdr": fdr,
                         "leading_edge": leading})
