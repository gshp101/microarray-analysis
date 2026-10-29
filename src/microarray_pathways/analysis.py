"""End-to-end analysis of one loaded dataset: genes -> metabolic pathways."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import de as de_mod
from . import enrichment, plots, report
from .geo import ExpressionDataset


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
) -> dict[str, pd.DataFrame]:
    """Run differential expression, ORA and GSEA; write tables, figures and report.html to out_dir/<accession>."""
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

    counts = {src: sum(k.startswith(f"{src}:") for k in gene_sets) for src in ("KEGG", "Reactome")}
    report.write_report(out, label, case, control, de, gsea, ora, ds.samples,
                        source=source or ds.accession, notes=ds.notes, fdr=fdr, min_abs_lfc=min_abs_lfc,
                        gene_set_counts={k: v for k, v in counts.items() if v})
    return {"genes": de, "ora": ora, "gsea": gsea}
