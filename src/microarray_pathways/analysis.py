"""End-to-end analysis of one loaded dataset: genes -> metabolic pathways."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import de as de_mod
from . import depmap, enrichment, plots, proximity, report, treatments
from .geo import DEFAULT_CACHE, ExpressionDataset


def analyze(
    ds: ExpressionDataset,
    case: str,
    control: str,
    gene_sets: dict[str, set[str]],
    out_dir: str | Path,
    label: str = "",
    fdr: float = 0.05,
    min_abs_lfc: float = 1.0,
    permutations: int = 1000,
    source: str = "",
    drugs: bool = True,
    cache_dir: Path = DEFAULT_CACHE,
    lineage: str | None = None,
) -> dict[str, pd.DataFrame]:
    """Run differential expression, ORA and GSEA, and (if `drugs`) the treatment-hypothesis steps;
    write tables, figures and report.html to out_dir/<accession>.

    `lineage` is the tissue of matching cancer cell lines (CCLE naming, e.g. "LUNG") for the
    DepMap dependency check; leave it out for non-cancer data."""
    out = Path(out_dir) / ds.accession
    out.mkdir(parents=True, exist_ok=True)
    label = label or f"{ds.accession}: {case} vs {control}"

    de = de_mod.run(ds, case, control)
    de.to_csv(out / "genes.tsv", sep="\t", float_format="%.4g")
    plots.volcano(de, out / "volcano.png", label, fdr=fdr, min_abs_lfc=min_abs_lfc)

    ora = enrichment.ora_up_down(de, gene_sets, fdr=fdr, min_abs_lfc=min_abs_lfc)
    ora.to_csv(out / "pathways_ora.tsv", sep="\t", index=False, float_format="%.4g")

    gsea = enrichment.gsea(de, gene_sets, permutations=permutations)
    gsea.to_csv(out / "pathways_gsea.tsv", sep="\t", index=False, float_format="%.4g")
    plots.top_pathways(gsea, out / "pathways_gsea.png", f"{label}\ntop metabolic pathways (GSEA)")

    n_sig = int(((de["fdr"] < fdr) & (de["log2fc"].abs() >= min_abs_lfc)).sum())
    summary = [
        label,
        f"samples: {case}={de.attrs['n_case']}, {control}={de.attrs['n_control']}",
        f"genes tested: {len(de)}; significant (FDR<{fdr}, |log2FC|>={min_abs_lfc}): {n_sig}",
        f"metabolic pathways tested (GSEA): {len(gsea)}; FDR<0.05: {int((gsea['fdr'] < 0.05).sum())}",
        f"ORA pathways with FDR<0.05: {int((ora['fdr'] < 0.05).sum())}",
    ]
    (out / "summary.txt").write_text("\n".join(summary) + "\n")

    drug_res, drug_error = None, ""
    if drugs:
        try:
            drug_res = treatments.treatment_hypotheses(de, gsea, gene_sets, out, cache_dir=cache_dir, fdr=fdr)
            summary.append(f"drugs acting on pathway genes (KEGG): {len(drug_res['targets'])}; "
                           f"opposing the change: {int(drug_res['targets']['opposes_change'].sum())}")
        except Exception as e:  # network or database problems should not lose the pathway report
            drug_error = f"{type(e).__name__}: {e}"
            summary.append(f"treatment hypotheses skipped: {drug_error}")
        if drug_res is not None:
            try:
                drug_res.update(proximity.proximity_hypotheses(
                    de, gsea, out, genes=drug_res["genes"], targets=drug_res["targets"],
                    reversal=drug_res["reversal"], cache_dir=cache_dir, fdr=fdr))
                prox = drug_res["proximity"]
                summary.append(f"drugs scored by network proximity: {len(prox)}; "
                               f"z < -2: {int((prox['z'] < -2).sum())}")
            except Exception as e:
                drug_res["proximity_error"] = f"{type(e).__name__}: {e}"
                summary.append(f"network proximity skipped: {drug_res['proximity_error']}")
            try:
                drug_res.update(depmap.dependency_hypotheses(drug_res["genes"], lineage, out,
                                                             targets=drug_res["targets"], cache_dir=cache_dir))
                dm, info = drug_res["depmap"], drug_res["depmap_info"]
                if "skipped" in info:
                    summary.append(f"DepMap dependencies skipped: {info['skipped']}")
                else:
                    summary.append(f"DepMap {info['lineage']} lines: {info['lineage_lines']}; pathway genes "
                                   f"they need more than other lines (FDR<0.1): {int((dm['fdr'] < 0.1).sum())}; "
                                   f"up and needed: {int(dm['up_and_needed'].sum())}")
            except Exception as e:
                drug_res["depmap_info"] = {"skipped": f"{type(e).__name__}: {e}"}
                summary.append(f"DepMap dependencies skipped: {drug_res['depmap_info']['skipped']}")
        (out / "summary.txt").write_text("\n".join(summary) + "\n")

    counts = {src: sum(k.startswith(f"{src}:") for k in gene_sets) for src in ("KEGG", "Reactome")}
    report.write_report(out, label, case, control, de, gsea, ora, ds.samples,
                        source=source or ds.accession, notes=ds.notes, fdr=fdr, min_abs_lfc=min_abs_lfc,
                        gene_set_counts={k: v for k, v in counts.items() if v},
                        drugs=drug_res, drug_error=drug_error)
    return {"genes": de, "ora": ora, "gsea": gsea, **({"drugs": drug_res} if drug_res else {})}
