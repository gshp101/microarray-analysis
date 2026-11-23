"""The PathRx pipeline as separate steps, for the desktop app.

Each method of Session is one step of the app. It calls the pathrx package exactly as
the command line does (native downloads, KEGG, gseapy) and keeps its output for the
next step. Downloads go through pathrx.geo.download first so the app can show their
progress; pathrx then finds them in its cache. Outputs are written like `pathrx report`
writes them: <out_dir>/<accession>/genes.tsv, pathways_gsea.tsv, report.html, ...

Nothing here imports Qt, so the steps can be tested and scripted on their own.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .. import datasets, depmap, enrichment, geo, pathways, plots, proximity, report, treatments
from .. import de as de_mod
from ..pathways import KEGG_REST

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# progress(label, loaded_bytes, total_bytes); total 0 = unknown, loaded -1 = busy, no bytes
Progress = Callable[[str, int, int], None]


@dataclass
class StepResult:
    """What a step shows: headline numbers, notes, figures and tables."""
    stats: list[tuple[str, str, str]] = field(default_factory=list)  # (value, label, tone: ""/up/down)
    notes: list[str] = field(default_factory=list)  # warnings shown above the output
    text: str = ""  # rich text (HTML) paragraph
    figures: list[tuple[str, Path]] = field(default_factory=list)  # (caption, png)
    tables: list[tuple[str, pd.DataFrame]] = field(default_factory=list)  # (tab title, table)
    data: dict = field(default_factory=dict)  # anything else the page needs


def _n(x) -> str:
    return f"{x:,}" if isinstance(x, (int, np.integer)) else str(x)


class Session:
    def __init__(self, cache_dir: Path = geo.DEFAULT_CACHE, out_dir: Path | str = "results",
                 progress: Progress | None = None):
        self.cache_dir = Path(cache_dir)
        self.out_root = Path(out_dir)
        self.progress = progress or (lambda *_: None)
        self.state: dict = {}

    # ------------------------------------------------------------ helpers

    @property
    def out(self) -> Path:
        out = self.out_root / self.state.get("accession", "pathrx")
        out.mkdir(parents=True, exist_ok=True)
        return out

    def busy(self, label: str) -> None:
        self.progress(label, -1, 0)

    def download(self, url: str, label: str = "") -> Path:
        label = label or url.rsplit("/", 1)[-1]
        if geo.cache_file(url, self.cache_dir).exists():
            return geo.cache_file(url, self.cache_dir)
        self.progress(label, 0, 0)
        return geo.download(url, self.cache_dir, progress=lambda done, total: self.progress(label, done, total))

    def prefetch(self, *urls: str, label: str = "") -> None:
        for u in urls:
            self.download(u, label)

    def _drop_after(self, keep: tuple[str, ...]) -> None:
        for k in [k for k in self.state if k not in keep]:
            self.state.pop(k)

    # ------------------------------------------------------------ step 1: find a dataset

    @staticmethod
    def presets() -> list[datasets.DatasetSpec]:
        return list(datasets.DATASETS.values())

    def search_geo(self, query: str, retmax: int = 20) -> tuple[int, list[dict]]:
        """Human expression-by-array series in GEO DataSets matching `query`, newest first."""
        import requests
        term = (f'({query}) AND gse[ETYP] AND "expression profiling by array"[DataSet Type] '
                f'AND "Homo sapiens"[Organism]')
        self.busy("Searching GEO")
        s = requests.get(f"{EUTILS}esearch.fcgi", params={"db": "gds", "retmode": "json", "retmax": retmax,
                                                          "term": term}, timeout=60)
        s.raise_for_status()
        res = s.json()["esearchresult"]
        ids = res.get("idlist", [])
        if not ids:
            return 0, []
        r = requests.get(f"{EUTILS}esummary.fcgi", params={"db": "gds", "retmode": "json", "id": ",".join(ids)},
                         timeout=60)
        r.raise_for_status()
        summ = r.json()["result"]
        hits = []
        for uid in summ.get("uids", []):
            d = summ[uid]
            hits.append({"accession": d.get("accession", ""), "title": d.get("title", ""),
                         "summary": d.get("summary", ""), "samples": int(d.get("n_samples") or 0),
                         "platforms": ", ".join("GPL" + g for g in str(d.get("gpl", "")).split(";") if g),
                         "date": d.get("pdat", "")})
        return int(res.get("count", len(hits))), hits

    # ------------------------------------------------------------ step 2: load the series

    def list_platforms(self, gse: str) -> list[dict]:
        """Series matrix files of a GSE, one per platform."""
        self.busy(f"Listing {gse} files")
        out = []
        for url in geo.series_matrix_urls(gse.strip().upper(), self.cache_dir):
            m = re.search(r"-(GPL\d+)_", url)
            out.append({"url": url, "file": url.rsplit("/", 1)[-1], "platform": m.group(1) if m else ""})
        return out

    def load_series(self, gse: str, url: str) -> StepResult:
        gse = gse.strip().upper()
        self.download(url, f"{gse} series matrix")
        self.busy("Reading the series matrix")
        expr, samples, gpl = geo.parse_series_matrix(geo.fetch(url, self.cache_dir).decode("utf-8", "replace"))
        annot = f"{geo.GEO_FTP}/platforms/{geo._prefix_dir(gpl)}/{gpl}/annot/{gpl}.annot.gz"
        try:
            self.download(annot, f"{gpl} annotation")
            sym_source = annot
        except Exception as e:  # some platforms have no curated annot file
            status = getattr(getattr(e, "response", None), "status_code", None)
            if status not in (403, 404):
                raise
            sym_source = f"{geo.GEO_ACC}?acc={gpl}&targ=self&form=text&view=data"
            self.download(sym_source, f"{gpl} platform table")
        self.busy("Mapping probes to genes")
        symbols = geo.parse_platform_table(geo.fetch(sym_source, self.cache_dir).decode("utf-8", "replace"))
        self.state = dict(accession=gse, url=url, raw_expr=expr, raw_samples=samples, platform=gpl,
                          symbols=symbols)
        mapped = symbols.reindex(expr.index).fillna("")
        single = int(((mapped != "") & ~mapped.str.contains("///")).sum())
        q = np.nanquantile(expr.to_numpy(dtype=float), [0, 0.25, 0.5, 0.75, 0.99, 1.0])
        log2 = geo.needs_log2(expr)
        meta = samples[[c for c in samples.columns if c.startswith(geo.META_PREFIXES)]].fillna("")
        plots.sample_distributions(expr, None, f"{gse}: values as deposited in GEO",
                                   self.out / "distribution_raw.png")
        self.state["sample_meta"] = meta
        qs = " / ".join(f"{x:.4g}" for x in q)
        return StepResult(
            stats=[(gse, "series", ""), (gpl, "platform", ""), (_n(expr.shape[1]), "samples", ""),
                   (_n(expr.shape[0]), "probes", ""), (_n(single), "probes with one gene symbol", "")],
            text=(f"Values min / 25% / median / 75% / 99% / max: <code>{qs}</code>. "
                  + ("They look linear, so the next step applies log2 (the GEO2R rule)." if log2
                     else "They already look log2-scaled.")
                  + f"<br>Gene symbols from <code>{sym_source}</code>"),
            figures=[("Distribution of values per sample", self.out / "distribution_raw.png")],
            tables=[("Sample descriptions", meta.rename_axis("sample").reset_index())],
            data={"fields": geo.metadata_fields(samples), "samples": int(expr.shape[1])},
        )

    # ------------------------------------------------------------ step 3: groups

    @staticmethod
    def _rules(case_name: str, case_rx: str, control_name: str, control_rx: str) -> dict[str, str]:
        case_name, control_name = case_name.strip() or "case", control_name.strip() or "control"
        if case_name == control_name:
            raise ValueError("the two groups need different names")
        return {case_name: case_rx, control_name: control_rx}

    def preview_groups(self, case_name, case_rx, control_name, control_rx, sample_filter="", exclude=""):
        """Group of each sample under the rules; fast enough to run on every keystroke."""
        rules = self._rules(case_name, case_rx, control_name, control_rx)
        return geo.group_status(self.state["raw_samples"], rules, sample_filter or None, exclude or None)

    def apply_groups(self, case_name, case_rx, control_name, control_rx, sample_filter="", exclude="") -> StepResult:
        rules = self._rules(case_name, case_rx, control_name, control_rx)
        if not case_rx or not control_rx:
            raise ValueError("write a pattern for both groups")
        self.busy("Building the gene table")
        s = self.state
        ds = geo.build_dataset(s["accession"], s["raw_expr"], s["raw_samples"], s["platform"], s["symbols"], rules,
                               sample_filter=sample_filter or None, exclude=exclude or None, source=s["url"])
        case, control = list(rules)
        counts = ds.group_counts()
        n_case, n_control = int(counts.get(case, 0)), int(counts.get(control, 0))
        if n_case < 2 or n_control < 2:
            raise ValueError(f"need at least 2 samples per group, got {case}={n_case}, {control}={n_control}")
        self._drop_after(("accession", "url", "raw_expr", "raw_samples", "platform", "symbols", "sample_meta"))
        s.update(ds=ds, case=case, control=control, rules=rules)
        plots.sample_distributions(ds.expression, ds.samples["group"],
                                   f"{ds.accession}: log2 expression per sample, after probe-to-gene collapse",
                                   self.out / "distribution_log2.png")
        return StepResult(
            stats=[(_n(n_case), case, "up"), (_n(n_control), control, "down"),
                   (_n(int(ds.samples["group"].isna().sum())), "left out", ""), (_n(ds.expression.shape[0]), "genes", "")],
            text="<br>".join(map(_esc, ds.notes)),
            figures=[("log2 expression per sample by group", self.out / "distribution_log2.png")],
            tables=[("Gene × sample table (log2)", ds.expression.round(3).rename_axis("gene").reset_index())],
        )

    # ------------------------------------------------------------ step 4: differential expression

    def run_de(self, fdr: float = 0.05, min_abs_lfc: float = 1.0) -> StepResult:
        s = self.state
        ds, case, control = s["ds"], s["case"], s["control"]
        self.busy("Testing every gene")
        de = de_mod.run(ds, case, control)
        self._drop_after(("accession", "url", "raw_expr", "raw_samples", "platform", "symbols", "sample_meta",
                          "ds", "case", "control", "rules"))
        s.update(de=de, fdr=fdr, min_abs_lfc=min_abs_lfc)
        de.to_csv(self.out / "genes.tsv", sep="\t", float_format="%.4g")
        plots.volcano(de, self.out / "volcano.png", f"{ds.accession}: {case} vs {control}", fdr=fdr,
                      min_abs_lfc=min_abs_lfc)
        sig = de_mod.significant(de, fdr, min_abs_lfc)
        return StepResult(
            stats=[(_n(len(de)), "genes tested", ""), (_n(int((sig["log2fc"] > 0).sum())), "up", "up"),
                   (_n(int((sig["log2fc"] < 0).sum())), "down", "down"),
                   (_n(int((de["fdr"] < fdr).sum())), f"FDR < {fdr:g}, any fold change", "")],
            figures=[("Volcano plot", self.out / "volcano.png")],
            tables=[("Genes", de.rename_axis("gene").reset_index())],
        )

    # ------------------------------------------------------------ step 5: gene sets

    def load_gene_sets(self, scope: str = "metabolic", sources: tuple[str, ...] = ("kegg", "reactome")) -> StepResult:
        if not sources:
            raise ValueError("pick at least one pathway database")
        if "reactome" in sources:
            self.prefetch(f"{pathways.REACTOME_DL}/ReactomePathways.gmt.zip",
                          f"{pathways.REACTOME_DL}/ReactomePathwaysRelation.txt", label="Reactome")
        if "kegg" in sources:
            urls = [f"{KEGG_REST}/list/pathway/hsa", f"{KEGG_REST}/link/hsa/pathway", f"{KEGG_REST}/list/hsa"]
            if scope != "metabolic":
                urls.append(f"{KEGG_REST}/get/br:br08901")
            self.prefetch(*urls, label="KEGG")
        self.busy("Reading gene sets")
        raw = pathways.pathway_gene_sets(tuple(sources), scope, self.cache_dir, max_machinery=None)
        sets = {k: v for k, v in raw.items() if pathways.machinery_fraction(v) <= 0.5}
        notes = []
        if len(raw) > len(sets):
            notes.append(f"{len(raw) - len(sets)} sets dropped: over half ribosome, proteasome, nuclear pore "
                         "or translation genes")
        measured = set(self.state["de"].index)
        df = pd.DataFrame([{"pathway": k, "source": k.split(":", 1)[0], "genes": len(v),
                            "measured": len(v & measured), "scored": 10 <= len(v & measured) <= 500}
                           for k, v in sets.items()],
                          columns=["pathway", "source", "genes", "measured", "scored"]).sort_values(["source", "pathway"])
        self._drop_after(("accession", "url", "raw_expr", "raw_samples", "platform", "symbols", "sample_meta",
                          "ds", "case", "control", "rules", "de", "fdr", "min_abs_lfc"))
        self.state.update(gene_sets=sets, scope=scope, kind={"all": "metabolic and signaling"}.get(scope, scope),
                          sources=sources)
        counts = df["source"].value_counts().to_dict()
        return StepResult(
            stats=[(_n(counts.get("Reactome", 0)), "Reactome sets", ""), (_n(counts.get("KEGG", 0)), "KEGG sets", ""),
                   (_n(int(df["scored"].sum())), "with 10 to 500 measured genes (scored)", "")],
            notes=notes, tables=[("Gene sets", df)],
        )

    # ------------------------------------------------------------ step 6: score pathways

    def run_pathways(self, permutations: int = 1000) -> StepResult:
        s = self.state
        de, sets, fdr, lfc, kind = s["de"], s["gene_sets"], s["fdr"], s["min_abs_lfc"], s["kind"]
        self.busy("Over-representation analysis")
        ora = enrichment.ora_up_down(de, sets, fdr=fdr, min_abs_lfc=lfc)
        engine = "gseapy" if _has_gseapy() else "numpy"
        self.busy(f"GSEA, {permutations} permutations ({engine})")
        gsea = enrichment.gsea(de, sets, permutations=permutations, engine=engine)
        self._drop_after(("accession", "url", "raw_expr", "raw_samples", "platform", "symbols", "sample_meta",
                          "ds", "case", "control", "rules", "de", "fdr", "min_abs_lfc", "gene_sets", "scope",
                          "kind", "sources"))
        s.update(ora=ora, gsea=gsea, permutations=permutations, gsea_engine=engine)
        out = self.out
        ora.to_csv(out / "pathways_ora.tsv", sep="\t", index=False, float_format="%.4g")
        gsea.to_csv(out / "pathways_gsea.tsv", sep="\t", index=False, float_format="%.4g")
        label = f"{s['ds'].accession}: {s['case']} vs {s['control']}"
        plots.top_pathways(gsea, out / "pathways_gsea.png", f"{label}\ntop {kind} pathways (GSEA)")
        sig = int((gsea["fdr"] < 0.05).sum())
        up = int(((gsea["fdr"] < 0.05) & (gsea["nes"] > 0)).sum())
        return StepResult(
            stats=[(_n(len(gsea)), "pathways scored", ""), (_n(sig), "GSEA FDR < 0.05", ""), (_n(up), "of those up", "up"),
                   (_n(sig - up), "down", "down"), (_n(int((ora["fdr"] < 0.05).sum())), "ORA FDR < 0.05", "")],
            text=report.findings(gsea, s["case"], s["control"], kind=kind),
            figures=[("Top pathways by GSEA", out / "pathways_gsea.png")],
            tables=[("GSEA", gsea), ("ORA (up and down genes)", ora)],
            data={"engine": engine},
        )

    # ------------------------------------------------------------ step 7: drug leads

    def run_drugs(self) -> StepResult:
        s = self.state
        de, gsea, sets, fdr = s["de"], s["gsea"], s["gene_sets"], s["fdr"]
        genes = treatments.pathway_genes(gsea, de, fdr=fdr)
        self.prefetch(f"{KEGG_REST}/list/hsa", f"{KEGG_REST}/link/drug/hsa", f"{KEGG_REST}/list/drug",
                      label="KEGG DRUG")
        gene_drugs = treatments.kegg_gene_drugs(self.cache_dir)
        wanted = sorted({d for g in genes.index for d in gene_drugs.get(g, ())})
        entries: dict[str, dict] = {}
        for i in range(0, len(wanted), 50):  # in chunks so progress shows; cached in one JSON file
            self.progress(f"KEGG DRUG entries {i + 1}-{min(i + 50, len(wanted))} of {len(wanted)}", i, len(wanted))
            entries.update(treatments.kegg_drug_entries(set(wanted[i:i + 50]), self.cache_dir))
        targets = treatments.target_candidates(genes, gene_drugs, entries)
        self.prefetch(treatments.ENRICHR_LIBRARY.format(treatments.L1000_LIBRARY), label="LINCS L1000 signatures")
        self.busy("Scoring L1000 signatures")
        metabolic = set().union(*sets.values()) if sets else set()
        reversal = treatments.reversal_candidates(de, treatments.l1000_signatures(self.cache_dir), metabolic, genes,
                                                  treatments.kegg_clinical_names(self.cache_dir))
        rev_z = {treatments._drug_key(c): z for c, z in zip(reversal["compound"], reversal["z_metabolic"])}
        targets["l1000_z_metabolic"] = [rev_z.get(treatments._drug_key(n), np.nan) for n in targets["drug"]]
        opposing = set(targets.loc[targets["opposes_change"].astype(bool), "drug"].map(treatments._drug_key))
        reversal["known_target_drug"] = reversal["compound"].map(treatments._drug_key).isin(opposing)
        out = self.out
        gtab = genes.assign(pathways=genes["pathways"].map("; ".join))
        gtab.to_csv(out / "drug_target_genes.tsv", sep="\t", float_format="%.4g", index_label="gene")
        targets.to_csv(out / "drugs_targets.tsv", sep="\t", index=False, float_format="%.4g")
        reversal.to_csv(out / "drugs_reversal.tsv", sep="\t", index=False, float_format="%.4g")
        s["drugs"] = {"genes": genes, "targets": targets, "reversal": reversal}
        n_opp = int(targets["opposes_change"].astype(bool).sum()) if len(targets) else 0
        return StepResult(
            stats=[(_n(len(genes)), "pathway-driving genes", ""),
                   (_n(int((genes["direction"] == "up").sum()) if len(genes) else 0), "of them up", "up"),
                   (_n(len(targets)), "KEGG drugs on them", ""), (_n(n_opp), "acting against the change", ""),
                   (_n(len(reversal)), "L1000 compounds scored", "")],
            text=("Signature reversal: negative z means the compound pushes these genes opposite to the disease. "
                  "Rank by percentile; z is a ranking score, not a test. These are research leads, not treatment "
                  "advice."),
            tables=[("Pathway-driving genes", gtab.rename_axis("gene").reset_index()),
                    ("Drugs on those genes (KEGG)", targets), ("Signature reversal (L1000)", reversal)],
        )

    # ------------------------------------------------------------ step 8: network proximity

    def run_proximity(self, reps: int = 1000) -> StepResult:
        drugs = self.state.get("drugs")
        if drugs is None:
            raise ValueError("run the drug leads step first")
        self.prefetch(proximity.INTERACTOME_URL, label="protein interaction network")
        self.prefetch(f"{KEGG_REST}/list/hsa", f"{KEGG_REST}/link/drug/hsa", f"{KEGG_REST}/list/drug", label="KEGG")
        self.busy(f"Network proximity, {reps} random reference sets")
        res = proximity.proximity_hypotheses(self.state["de"], self.state["gsea"], self.out, genes=drugs["genes"],
                                             targets=drugs["targets"], reversal=drugs["reversal"],
                                             cache_dir=self.cache_dir, fdr=self.state["fdr"], reps=reps)
        drugs.pop("proximity_error", None)
        drugs.update(res)
        prox, info = res["proximity"], res["proximity_info"]
        clin = prox[prox["clinical"].astype(bool)] if len(prox) else prox
        return StepResult(
            stats=[(_n(info["module_genes"]), f"of {info['pathway_genes']} pathway genes in the network", ""),
                   (_n(info["network_nodes"]), "proteins in the network", ""), (_n(len(prox)), "drugs scored", ""),
                   (_n(len(clin)), "clinical drugs", ""),
                   (_n(int((clin["z"] < -2).sum()) if len(clin) else 0), "clinical with z < −2", "")],
            tables=[("Drugs by network proximity", prox)],
        )

    # ------------------------------------------------------------ step 9: DepMap

    def _depmap_libs(self) -> None:
        libs = [depmap.DEPENDENCY_LIBRARY, *depmap.TISSUE_LIBRARIES]
        self.prefetch(*(treatments.ENRICHR_LIBRARY.format(lib) for lib in libs), label="DepMap via Enrichr")

    def lineages(self) -> list[tuple[str, int]]:
        """Tissues of the CCLE cell lines DepMap covers, with line counts."""
        self._depmap_libs()
        deps, tissues = depmap.load_dependencies(self.cache_dir), depmap.load_tissues(self.cache_dir)
        counts = pd.Series([tissues[c] for c in deps if c in tissues]).value_counts()
        return [(str(k), int(v)) for k, v in counts.items()]

    def run_depmap(self, lineage: str | None) -> StepResult:
        drugs = self.state.get("drugs")
        if drugs is None:
            raise ValueError("run the drug leads step first")
        lineage = (lineage or "").strip() or None
        if lineage:
            self._depmap_libs()
        self.busy("Checking cell line dependencies")
        res = depmap.dependency_hypotheses(drugs["genes"], lineage, self.out, targets=drugs["targets"],
                                           cache_dir=self.cache_dir)
        drugs.update(res)
        info, dm = res["depmap_info"], res["depmap"]
        if "skipped" in info:
            return StepResult(text=f"Skipped: {_esc(info['skipped'])}.", data={"skipped": True})
        return StepResult(
            stats=[(_n(info["lineage_lines"]), f"{info['lineage'].lower()} cell lines", ""),
                   (_n(info["lines"]), "cell lines with a known tissue", ""),
                   (_n(int((dm["fdr"] < 0.1).sum()) if len(dm) else 0), "genes needed more often (FDR < 0.1)", ""),
                   (_n(int(dm["up_and_needed"].sum()) if len(dm) else 0), "up and needed", "up")],
            tables=[("DepMap dependencies", dm)] if len(dm) else [],
            text="" if len(dm) else "No pathway-driving genes to check.",
        )

    # ------------------------------------------------------------ step 10: report

    def make_report(self) -> Path:
        s = self.state
        ds, case, control, sets = s["ds"], s["case"], s["control"], s["gene_sets"]
        self.busy("Writing the report")
        counts = {src: sum(k.startswith(f"{src}:") for k in sets) for src in ("KEGG", "Reactome")}
        notes = list(ds.notes) + [f"run in the PathRx desktop app; GSEA by {s.get('gsea_engine', 'gseapy')}, "
                                  f"{s.get('permutations', 1000)} permutations"]
        drugs = s.get("drugs")
        source = f"GEO {ds.accession}"
        spec = next((d for d in datasets.DATASETS.values() if d.accession == ds.accession), None)
        if spec:
            source += f" ({spec.description})"
        summary = [f"{ds.accession}: {case} vs {control}",
                   f"samples: {case}={s['de'].attrs['n_case']}, {control}={s['de'].attrs['n_control']}",
                   f"{s['kind']} pathways tested (GSEA): {len(s['gsea'])}; "
                   f"FDR<0.05: {int((s['gsea']['fdr'] < 0.05).sum())}"]
        (self.out / "summary.txt").write_text("\n".join(summary) + "\n")
        return report.write_report(self.out, f"{ds.accession}: {case} vs {control}", case, control, s["de"],
                                   s["gsea"], s["ora"], ds.samples, source=source, notes=notes, fdr=s["fdr"],
                                   min_abs_lfc=s["min_abs_lfc"],
                                   gene_set_counts={k: v for k, v in counts.items() if v}, drugs=drugs,
                                   drug_error="" if drugs else "the drug leads step was not run", kind=s["kind"])


def _has_gseapy() -> bool:
    try:
        import gseapy  # noqa: F401
        return True
    except ImportError:
        return False


def _esc(x) -> str:
    import html
    return html.escape(str(x))
