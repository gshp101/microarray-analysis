"""Which pathway-driving genes do cancer cell lines of the matching tissue depend on?

DepMap knocks out each gene with CRISPR in ~1,100 cancer cell lines and measures
whether the cells keep growing. A gene that went up in the tumors, drives a shifted
metabolic pathway, and is needed by cell lines of the same tissue more often than
by other cell lines is a better target lead than one that merely changed.

Data (the DepMap portal and figshare are not reachable from this project's
environment, so both come through Enrichr):
  dependencies  DepMap_CRISPR_GeneDependency_CellLines_2023: for each cell line, the
                genes it depends on, minus the common essential genes every line needs
                (ribosome, polymerase and the like), so what is left are selective
                dependencies.
  tissue        CCLE cell line names, which end in the tissue of origin
                (A549_LUNG), from the Enrichr CCLE expression and proteomics libraries.

For each gene: the share of tissue-matched lines that depend on it, the share of all
other lines, and a one-sided Fisher exact test that the tissue lines depend on it
more often.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from .geo import DEFAULT_CACHE, fetch
from .treatments import ENRICHR_LIBRARY

DEPENDENCY_LIBRARY = "DepMap_CRISPR_GeneDependency_CellLines_2023"
TISSUE_LIBRARIES = ("Cancer_Cell_Line_Encyclopedia", "CCLE_Proteomics_2020")


def norm_line(name: str) -> str:
    """'NCI-H1299', 'NCIH1299' and 'NCI H1299' -> 'NCIH1299'."""
    return re.sub(r"[^A-Z0-9]", "", name.upper())


def parse_gmt(text: str) -> dict[str, set[str]]:
    sets = {}
    for line in text.splitlines():
        row = line.split("\t")
        if len(row) >= 3:
            sets[row[0].strip()] = {g.split(",")[0] for g in row[2:] if g}
    return sets


def parse_tissues(terms) -> dict[str, str]:
    """CCLE terms ('A549 LUNG', '22RV1 PROSTATE TenPx24') -> {normalized line: tissue}."""
    out = {}
    for term in terms:
        parts = re.sub(r"\s+TenPx\d+$", "", term.strip()).split(" ", 1)
        if len(parts) == 2:
            out.setdefault(norm_line(parts[0]), parts[1].upper())
    return out


def load_dependencies(cache_dir: Path = DEFAULT_CACHE) -> dict[str, set[str]]:
    sets = parse_gmt(fetch(ENRICHR_LIBRARY.format(DEPENDENCY_LIBRARY), cache_dir).decode())
    merged: dict[str, set[str]] = {}
    for name, genes in sets.items():
        merged.setdefault(norm_line(name), set()).update(genes)
    return merged


def load_tissues(cache_dir: Path = DEFAULT_CACHE) -> dict[str, str]:
    out: dict[str, str] = {}
    for lib in TISSUE_LIBRARIES:
        for k, v in parse_tissues(parse_gmt(fetch(ENRICHR_LIBRARY.format(lib), cache_dir).decode())).items():
            out.setdefault(k, v)
    return out


def lineage_dependencies(genes: pd.DataFrame, deps: dict[str, set[str]], tissues: dict[str, str],
                         lineage: str, drugs: dict[str, list[str]] | None = None) -> pd.DataFrame:
    """Dependency of tissue-matched vs other cell lines on each gene in `genes`.

    `genes` comes from treatments.pathway_genes (log2fc, direction, pathways). Only
    cell lines with a known tissue are used. `drugs` maps gene -> KEGG drugs that act
    against its change, to list beside it.
    """
    lines = [c for c in deps if c in tissues]
    in_lin = [c for c in lines if tissues[c] == lineage.upper()]
    other = [c for c in lines if tissues[c] != lineage.upper()]
    cols = ["gene", "direction", "log2fc", "lineage_dependent", "lineage_lines", "frac_lineage",
            "frac_other", "odds_ratio", "pvalue", "fdr", "up_and_needed", "drugs", "pathways"]
    if not in_lin or genes.empty:
        return pd.DataFrame(columns=cols)
    rows = []
    for g in genes.index:
        a = sum(g in deps[c] for c in in_lin)
        b = sum(g in deps[c] for c in other)
        odds, p = fisher_exact([[a, len(in_lin) - a], [b, len(other) - b]], alternative="greater")
        pw = genes.at[g, "pathways"]
        rows.append({"gene": g, "direction": genes.at[g, "direction"], "log2fc": genes.at[g, "log2fc"],
                     "lineage_dependent": a, "lineage_lines": len(in_lin), "frac_lineage": a / len(in_lin),
                     "frac_other": b / len(other) if other else np.nan, "odds_ratio": odds, "pvalue": p,
                     "drugs": ", ".join((drugs or {}).get(g, [])),
                     "pathways": "; ".join(pw) if isinstance(pw, (list, tuple, set)) else str(pw)})
    res = pd.DataFrame(rows)
    res["fdr"] = multipletests(res["pvalue"], method="fdr_bh")[1]
    # went up in the disease, and at least 1 in 10 tissue-matched lines needs it, more often than other lines
    res["up_and_needed"] = ((res["direction"] == "up") & (res["frac_lineage"] >= 0.1)
                            & (res["frac_lineage"] > res["frac_other"]))
    return res[cols].sort_values(["pvalue", "frac_lineage"], ascending=[True, False]).reset_index(drop=True)


def opposing_drugs_by_gene(targets: pd.DataFrame | None) -> dict[str, list[str]]:
    """Gene -> KEGG drugs that act against its change, from the drug-target table."""
    out: dict[str, list[str]] = {}
    if targets is None or targets.empty:
        return out
    for r in targets[targets["opposes_change"].astype(bool)].itertuples():
        for g in re.findall(r"([A-Za-z0-9\-_.]+) \(", str(r.targets)):
            out.setdefault(g, []).append(r.drug)
    return out


def dependency_hypotheses(genes: pd.DataFrame, lineage: str | None, out_dir: str | Path,
                          targets: pd.DataFrame | None = None, cache_dir: Path = DEFAULT_CACHE) -> dict:
    """Run the lineage dependency check and write depmap_dependencies.tsv.

    Returns {"depmap": table, "depmap_info": {...}}; with no lineage (non-cancer data, or
    a tissue DepMap does not cover) the table is empty and info says why.
    """
    if not lineage:
        return {"depmap": pd.DataFrame(), "depmap_info": {"skipped": "no matching cancer cell line tissue"}}
    deps, tissues = load_dependencies(cache_dir), load_tissues(cache_dir)
    res = lineage_dependencies(genes, deps, tissues, lineage, opposing_drugs_by_gene(targets))
    lines = [c for c in deps if c in tissues]
    info = {"lineage": lineage.upper(), "lineage_lines": sum(tissues[c] == lineage.upper() for c in lines),
            "lines": len(lines), "all_lines": len(deps)}
    if info["lineage_lines"] == 0:
        info["skipped"] = f"no DepMap cell lines with tissue {lineage.upper()}"
    res.to_csv(Path(out_dir) / "depmap_dependencies.tsv", sep="\t", index=False, float_format="%.4g")
    return {"depmap": res, "depmap_info": info}
