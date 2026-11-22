"""Browser glue for PathRx: runs the pipeline step by step inside Pyodide.

The pathrx package does the analysis exactly as the command line does. This module
only does what the browser needs on top of it:

  * downloads: pathrx.geo.fetch() reads from an on-disk cache, so each step first
    downloads its URLs with the browser's fetch (via `js_fetch`, given by the worker)
    into that cache, then calls the normal pathrx functions;
  * KEGG: rest.kegg.jp sends no CORS headers, so browsers block it. KEGG URLs are
    rewritten to a proxy prefix the user configures; with none, KEGG steps are skipped
    and the result says so;
  * results go back to the page as JSON (table previews, PNG plots as base64, notes).

Each public coroutine is one step and keeps its output in STATE for the next step.
"""

from __future__ import annotations

import base64
import io
import json
import re
import sys
import types
from collections import defaultdict
from pathlib import Path

# requests is imported at module level by pathrx; in the browser every download goes
# through js_fetch instead, so a stub is enough.
if "requests" not in sys.modules:
    _req = types.ModuleType("requests")

    class HTTPError(Exception):
        pass

    def _no_requests(*_a, **_k):
        raise RuntimeError("network calls in the browser go through pathrx_web.download")

    _req.HTTPError, _req.get = HTTPError, _no_requests
    sys.modules["requests"] = _req

import numpy as np
import pandas as pd

from pathrx import de as de_mod
from pathrx import depmap, enrichment, geo, pathways, plots, proximity, report, treatments
from pathrx.pathways import KEGG_REST

CACHE = Path("/cache")
OUT = Path("/out")
STATE: dict = {}
SETTINGS = {"kegg_proxy": ""}
js_fetch = None  # set by the worker: async (url, label) -> Uint8Array, raises on HTTP errors
MAX_ROWS = 1000


class DownloadError(Exception):
    def __init__(self, url: str, status: int, message: str = ""):
        super().__init__(message or f"HTTP {status} for {url}")
        self.url, self.status = url, status


# ------------------------------------------------------------ downloads

def cache_path(url: str) -> Path:
    """Same file name pathrx.geo.fetch() uses, so it finds what we download."""
    return CACHE / re.sub(r"[^A-Za-z0-9._-]", "_", url.split("://", 1)[-1])


def kegg_available() -> bool:
    return bool(SETTINGS["kegg_proxy"])


def _real_url(url: str) -> str:
    if url.startswith(KEGG_REST):
        proxy = SETTINGS["kegg_proxy"]
        if not proxy:
            raise DownloadError(url, 0, "KEGG blocks browser requests (no CORS); set a KEGG proxy to use it")
        return proxy.rstrip("/") + "/" + url[len(KEGG_REST):].lstrip("/")
    return url


async def download(url: str, label: str = "", cache: bool = True) -> bytes:
    """Bytes of `url` (raw, still gzipped if it is). Cached under pathrx's cache name."""
    path = cache_path(url)
    if cache and path.exists():
        return path.read_bytes()
    try:
        arr = await js_fetch(_real_url(url), label or url.rsplit("/", 1)[-1])
    except Exception as e:  # JS errors carry "HTTP <status>" in the message
        m = re.search(r"HTTP (\d+)", str(e))
        raise DownloadError(url, int(m.group(1)) if m else 0, str(e)) from None
    data = arr.to_bytes()
    if cache:
        CACHE.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return data


async def prefetch(*urls: str, label: str = "") -> None:
    for u in urls:
        await download(u, label)


def text_of(data: bytes) -> str:
    import gzip
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8", "replace")


# ------------------------------------------------------------ result helpers

def _cell(v):
    if isinstance(v, (set, frozenset, list, tuple)):
        return ", ".join(map(str, sorted(v) if isinstance(v, (set, frozenset)) else v))
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return None if not np.isfinite(v) else float(v)
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    return str(v) if not isinstance(v, (int, str)) else v


def table(df: pd.DataFrame, name: str, index: str | None = None, max_rows: int = MAX_ROWS) -> dict:
    """Preview of a table for the page, and keep the full table for download as name.tsv."""
    if index is not None:
        df = df.rename_axis(index).reset_index()
    STATE.setdefault("tables", {})[name] = df
    head = df.head(max_rows)
    return {"name": name, "columns": [str(c) for c in df.columns], "total": int(len(df)),
            "rows": [[_cell(v) for v in row] for row in head.itertuples(index=False)]}


def png(path: Path) -> str:
    return base64.b64encode(Path(path).read_bytes()).decode()


def fig_png(fig) -> str:
    import matplotlib.pyplot as plt
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def table_tsv(name: str) -> str:
    df = STATE["tables"][name]
    return df.to_csv(sep="\t", index=False, float_format="%.4g")


def configure(settings_json: str) -> str:
    SETTINGS.update(json.loads(settings_json))
    return json.dumps({"kegg": kegg_available()})


def dumps(obj) -> str:
    return json.dumps(obj, default=_cell)


# ------------------------------------------------------------ step 1: load a GEO series

META_PREFIXES = ("title", "source_name", "characteristics", "description")


async def list_platforms(gse: str) -> str:
    """Series matrix files of a GSE, one per platform."""
    gse = gse.strip().upper()
    base = f"{geo.GEO_FTP}/series/{geo._prefix_dir(gse)}/{gse}/matrix/"
    listing = text_of(await download(base, "file list", cache=False))
    names = sorted(set(re.findall(r'href="([^"]+_series_matrix\.txt\.gz)"', listing)))
    if not names:
        raise FileNotFoundError(f"no series matrix files under {base}")
    out = []
    for n in names:
        m = re.search(r"-(GPL\d+)_", n)
        out.append({"url": base + n, "file": n, "platform": m.group(1) if m else ""})
    return dumps({"accession": gse, "files": out})


async def load_series(gse: str, url: str) -> str:
    """Download and parse the series matrix and its platform annotation."""
    gse = gse.strip().upper()
    raw = await download(url, f"{gse} series matrix")
    expr, samples, gpl = geo.parse_series_matrix(text_of(raw))
    annot = f"{geo.GEO_FTP}/platforms/{geo._prefix_dir(gpl)}/{gpl}/annot/{gpl}.annot.gz"
    try:
        symbols = geo.parse_platform_table(text_of(await download(annot, f"{gpl} annotation")))
        sym_source = annot
    except DownloadError as e:
        if e.status not in (403, 404):
            raise
        full = f"{geo.GEO_ACC}?acc={gpl}&targ=self&form=text&view=data"
        symbols = geo.parse_platform_table(text_of(await download(full, f"{gpl} platform table")))
        sym_source = full
    STATE.clear()
    STATE.update(accession=gse, url=url, raw_expr=expr, raw_samples=samples, platform=gpl, symbols=symbols)

    mapped = symbols.reindex(expr.index).fillna("")
    single = int(((mapped != "") & ~mapped.str.contains("///")).sum())
    vals = expr.to_numpy(dtype=float)
    q = np.nanquantile(vals, [0, 0.25, 0.5, 0.75, 0.99, 1.0])
    meta_cols = [c for c in samples.columns if c.startswith(META_PREFIXES)]
    meta = samples[meta_cols].fillna("")
    # columns whose values repeat are the useful ones for defining groups
    fields = []
    for c in meta_cols:
        counts = meta[c].value_counts()
        if 1 < len(counts) <= 12 or (len(counts) == 1 and c.startswith("characteristics")):
            fields.append({"column": c, "values": [[str(k), int(v)] for k, v in counts.items()]})
    return dumps({
        "accession": gse, "platform": gpl, "samples": int(expr.shape[1]), "probes": int(expr.shape[0]),
        "mapped_probes": single, "symbol_source": sym_source,
        "quantiles": [float(x) for x in q], "needs_log2": geo.needs_log2(expr),
        "sample_table": table(meta, "samples", index="sample"),
        "fields": fields,
        "distribution": _distribution_plot(expr, None, f"{gse}: values as deposited in GEO"),
    })


def _distribution_plot(expr: pd.DataFrame, groups: pd.Series | None, title: str, max_samples: int = 120) -> str:
    import matplotlib.pyplot as plt
    cols = list(expr.columns)
    if len(cols) > max_samples:
        cols = cols[:: int(np.ceil(len(cols) / max_samples))]
    data = [expr[c].dropna().to_numpy() for c in cols]
    colors = []
    palette = {0: plots.UP, 1: plots.DOWN}
    names = list(dict.fromkeys(groups.dropna())) if groups is not None else []
    for c in cols:
        g = groups.get(c) if groups is not None else None
        colors.append(palette.get(names.index(g), plots.NEUTRAL) if g in names else plots.NEUTRAL)
    fig, ax = plt.subplots(figsize=(8, 3.2), dpi=120)
    plots._style(ax)
    bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.7)
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_edgecolor(plots.MUTED)
        patch.set_linewidth(0.4)
    for k in ("whiskers", "caps", "medians"):
        for line in bp[k]:
            line.set_color(plots.MUTED)
            line.set_linewidth(0.5)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(cols)} samples" + (f" (every {len(expr.columns) // len(cols)}th shown)"
                                            if len(cols) < len(expr.columns) else ""), color=plots.INK, fontsize=9)
    ax.set_ylabel("expression", color=plots.INK, fontsize=9)
    ax.set_title(title, loc="left", fontsize=10, color=plots.INK)
    if names:
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=palette[i], label=n) for i, n in enumerate(names[:2])],
                  frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    return fig_png(fig)


# ------------------------------------------------------------ step 2: groups

def _rules(case_name, case_rx, control_name, control_rx):
    case_name, control_name = case_name.strip() or "case", control_name.strip() or "control"
    if case_name == control_name:
        raise ValueError("the two groups need different names")
    for rx in (case_rx, control_rx):
        re.compile(rx)
    return {case_name: case_rx, control_name: control_rx}


def preview_groups(case_name: str, case_rx: str, control_name: str, control_rx: str,
                   sample_filter: str = "", exclude: str = "") -> str:
    """Group of each sample under the rules (fast; run while the user types)."""
    samples = STATE["raw_samples"]
    try:
        rules = _rules(case_name, case_rx, control_name, control_rx)
        text = geo.sample_text(samples)
        keep = pd.Series(True, index=samples.index)
        if sample_filter:
            keep &= text.str.contains(sample_filter, case=False, regex=True)
        if exclude:
            keep &= ~text.str.contains(exclude, case=False, regex=True)
        group = geo.assign_groups(samples, rules) if case_rx and control_rx else pd.Series(np.nan, index=samples.index)
        hits = pd.DataFrame({g: text.str.contains(rx, case=False, regex=True) if rx else False
                             for g, rx in rules.items()})
    except (re.error, ValueError) as e:
        return dumps({"error": str(e)})
    both = hits.sum(axis=1) > 1
    status = np.where(~keep, "filtered out", np.where(both, "matches both", group.fillna("unassigned")))
    counts = pd.Series(status).value_counts().to_dict()
    return dumps({"status": dict(zip(samples.index, status)), "counts": counts})


def apply_groups(case_name: str, case_rx: str, control_name: str, control_rx: str,
                 sample_filter: str = "", exclude: str = "") -> str:
    """Build the gene x sample dataset exactly as pathrx load does."""
    rules = _rules(case_name, case_rx, control_name, control_rx)
    ds = geo.build_dataset(STATE["accession"], STATE["raw_expr"], STATE["raw_samples"], STATE["platform"],
                           STATE["symbols"], rules, sample_filter=sample_filter or None,
                           exclude=exclude or None, source=STATE["url"])
    case, control = list(rules)
    counts = ds.group_counts()
    n_case, n_control = int(counts.get(case, 0)), int(counts.get(control, 0))
    if n_case < 2 or n_control < 2:
        raise ValueError(f"need at least 2 samples per group, got {case}={n_case}, {control}={n_control}")
    for k in [k for k in STATE if k not in ("accession", "url", "raw_expr", "raw_samples", "platform", "symbols")]:
        STATE.pop(k)
    STATE.update(ds=ds, case=case, control=control, rules=rules)
    e = ds.expression
    return dumps({
        "notes": ds.notes, "genes": int(e.shape[0]), "case": case, "control": control,
        "n_case": n_case, "n_control": n_control, "unassigned": int(ds.samples["group"].isna().sum()),
        "distribution": _distribution_plot(e, ds.samples["group"], f"{ds.accession}: log2 expression per sample, "
                                                                     "after probe-to-gene collapse"),
        "expression": table(e.round(3), "expression", index="gene", max_rows=200),
    })


# ------------------------------------------------------------ step 3: differential expression

def run_de(fdr: float = 0.05, min_abs_lfc: float = 1.0) -> str:
    ds, case, control = STATE["ds"], STATE["case"], STATE["control"]
    de = de_mod.run(ds, case, control)
    STATE.update(de=de, fdr=fdr, min_abs_lfc=min_abs_lfc)
    OUT.mkdir(parents=True, exist_ok=True)
    label = f"{ds.accession}: {case} vs {control}"
    plots.volcano(de, OUT / "volcano.png", label, fdr=fdr, min_abs_lfc=min_abs_lfc)
    sig = de_mod.significant(de, fdr, min_abs_lfc)
    return dumps({
        "tested": int(len(de)), "up": int((sig["log2fc"] > 0).sum()), "down": int((sig["log2fc"] < 0).sum()),
        "fdr_only": int((de["fdr"] < fdr).sum()), "n_case": de.attrs["n_case"], "n_control": de.attrs["n_control"],
        "volcano": png(OUT / "volcano.png"),
        "genes": table(de, "genes", index="gene"),
    })


# ------------------------------------------------------------ step 4: pathway gene sets

async def load_gene_sets(scope: str = "metabolic", use_kegg: bool = True) -> str:
    notes = []
    sets: dict[str, set[str]] = {}
    sources = []
    rx = f"{pathways.REACTOME_DL}/ReactomePathways.gmt.zip"
    rel = f"{pathways.REACTOME_DL}/ReactomePathwaysRelation.txt"
    await prefetch(rx, rel, label="Reactome")
    sets.update(pathways.reactome_sets(scope, CACHE))
    sources.append("reactome")
    if use_kegg and kegg_available():
        try:
            urls = [f"{KEGG_REST}/list/pathway/hsa", f"{KEGG_REST}/link/hsa/pathway", f"{KEGG_REST}/list/hsa"]
            if scope != "metabolic":
                urls.append(f"{KEGG_REST}/get/br:br08901")
            await prefetch(*urls, label="KEGG")
            sets.update(pathways.kegg_sets(scope, CACHE))
            sources.append("kegg")
        except DownloadError as e:
            notes.append(f"KEGG gene sets skipped: {e}")
    else:
        notes.append("KEGG gene sets skipped: KEGG blocks browser requests; set a KEGG proxy to include them. "
                     "Results use Reactome only, so they will differ from the command line, which uses both.")
    before = len(sets)
    sets = {k: v for k, v in sets.items() if pathways.machinery_fraction(v) <= 0.5}
    if before > len(sets):
        notes.append(f"{before - len(sets)} sets dropped: over half ribosome, proteasome, nuclear pore "
                     "or translation genes")
    measured = set(STATE["de"].index) if "de" in STATE else set()
    rows = [{"pathway": k, "source": k.split(":", 1)[0], "genes": len(v),
             "measured": len(v & measured) if measured else None,
             "scored": (10 <= len(v & measured) <= 500) if measured else None} for k, v in sets.items()]
    df = pd.DataFrame(rows).sort_values(["source", "pathway"])
    STATE.update(gene_sets=sets, scope=scope, kind={"all": "metabolic and signaling"}.get(scope, scope),
                 sources=sources)
    counts = df["source"].value_counts().to_dict()
    return dumps({"counts": counts, "total": len(sets), "notes": notes,
                  "scored": int(df["scored"].sum()) if measured else None,
                  "sets": table(df, "gene_sets")})


# ------------------------------------------------------------ step 5: pathway scoring

def run_pathways(permutations: int = 1000) -> str:
    de, sets, fdr, lfc = STATE["de"], STATE["gene_sets"], STATE["fdr"], STATE["min_abs_lfc"]
    ds, kind = STATE["ds"], STATE["kind"]
    ora = enrichment.ora_up_down(de, sets, fdr=fdr, min_abs_lfc=lfc)
    gsea = enrichment.gsea(de, sets, permutations=permutations, engine="numpy")
    STATE.update(ora=ora, gsea=gsea, permutations=permutations)
    label = f"{ds.accession}: {STATE['case']} vs {STATE['control']}"
    plots.top_pathways(gsea, OUT / "pathways_gsea.png", f"{label}\ntop {kind} pathways (GSEA)")
    return dumps({
        "gsea_tested": int(len(gsea)), "gsea_sig": int((gsea["fdr"] < 0.05).sum()),
        "gsea_up": int(((gsea["fdr"] < 0.05) & (gsea["nes"] > 0)).sum()),
        "ora_sig": int((ora["fdr"] < 0.05).sum()),
        "findings": report.findings(gsea, STATE["case"], STATE["control"], kind=kind),
        "chart": png(OUT / "pathways_gsea.png"),
        "gsea": table(gsea, "pathways_gsea"), "ora": table(ora, "pathways_ora"),
    })


# ------------------------------------------------------------ step 6: drug leads

TARGET_COLS = ["drug", "kegg_id", "action", "opposes_change", "marketed", "targets", "n_targets", "score",
               "pathways", "efficacy", "atc"]


async def _kegg_drug_entries(ids: set[str]) -> dict[str, dict]:
    cache = STATE.setdefault("kegg_drug_cache", {})
    todo = sorted(set(ids) - set(cache))
    for i in range(0, len(todo), 10):
        chunk = todo[i:i + 10]
        url = f"{KEGG_REST}/get/" + "+".join(f"dr:{d}" for d in chunk)
        text = text_of(await download(url, f"KEGG DRUG entries {i + 1}-{i + len(chunk)} of {len(todo)}"))
        for part in text.split("///"):
            if part.strip():
                d = treatments.parse_kegg_drug(part.strip() + "\n///")
                cache[d["id"]] = d
        for d in chunk:
            cache.setdefault(d, None)
    return {d: cache[d] for d in ids if cache.get(d)}


async def run_drugs() -> str:
    de, gsea, sets, fdr = STATE["de"], STATE["gsea"], STATE["gene_sets"], STATE["fdr"]
    genes = treatments.pathway_genes(gsea, de, fdr=fdr)
    notes = []
    targets = pd.DataFrame(columns=TARGET_COLS)
    clinical: set[str] = set()
    if kegg_available():
        await prefetch(f"{KEGG_REST}/list/hsa", f"{KEGG_REST}/link/drug/hsa", f"{KEGG_REST}/list/drug",
                       label="KEGG DRUG")
        gene_drugs = treatments.kegg_gene_drugs(CACHE)
        wanted = {d for g in genes.index for d in gene_drugs.get(g, ())}
        targets = treatments.target_candidates(genes, gene_drugs, await _kegg_drug_entries(wanted))
        if targets.empty:
            targets = pd.DataFrame(columns=TARGET_COLS)
        clinical = treatments.kegg_clinical_names(CACHE)
    else:
        notes.append("Drug-target route skipped: it needs KEGG DRUG, which blocks browser requests. "
                     "Set a KEGG proxy to run it. The L1000 route below still runs, but without KEGG it "
                     "cannot flag which compounds are clinical drugs.")
    await prefetch(treatments.ENRICHR_LIBRARY.format(treatments.L1000_LIBRARY), label="LINCS L1000 signatures")
    sigs = treatments.l1000_signatures(CACHE)
    metabolic = set().union(*sets.values()) if sets else set()
    reversal = treatments.reversal_candidates(de, sigs, metabolic, genes, clinical)
    rev_z = {treatments._drug_key(c): z for c, z in zip(reversal["compound"], reversal["z_metabolic"])}
    targets["l1000_z_metabolic"] = [rev_z.get(treatments._drug_key(n), np.nan) for n in targets["drug"]]
    opposing = set(targets.loc[targets["opposes_change"].astype(bool), "drug"].map(treatments._drug_key))
    reversal["known_target_drug"] = reversal["compound"].map(treatments._drug_key).isin(opposing)
    STATE["drugs"] = {"genes": genes, "targets": targets, "reversal": reversal}
    if not kegg_available():
        STATE["drug_note"] = notes[0]
    gtab = genes.assign(pathways=genes["pathways"].map("; ".join))
    return dumps({
        "notes": notes, "kegg": kegg_available(), "n_genes": int(len(genes)),
        "n_up": int((genes["direction"] == "up").sum()) if len(genes) else 0,
        "n_targets": int(len(targets)), "n_opposing": int(targets["opposes_change"].astype(bool).sum()),
        "n_compounds": int(len(reversal)), "n_signatures": len(sigs),
        "genes": table(gtab, "drug_target_genes", index="gene"),
        "targets": table(targets, "drugs_targets"),
        "reversal": table(reversal, "drugs_reversal"),
    })


# ------------------------------------------------------------ step 7: network proximity

async def run_proximity(reps: int = 1000) -> str:
    drugs = STATE.get("drugs")
    if drugs is None:
        raise ValueError("run the drug leads step first")
    if not kegg_available():
        drugs.pop("proximity", None)
        drugs["proximity_error"] = "needs KEGG (drug targets and gene ids), which blocks browser requests"
        return dumps({"skipped": "Network proximity needs KEGG DRUG targets and KEGG's gene id list, and KEGG "
                                 "blocks browser requests. Set a KEGG proxy to run this step."})
    await prefetch(proximity.INTERACTOME_URL, label="protein interaction network")
    await prefetch(f"{KEGG_REST}/list/hsa", f"{KEGG_REST}/link/drug/hsa", f"{KEGG_REST}/list/drug", label="KEGG")
    res = proximity.proximity_hypotheses(STATE["de"], STATE["gsea"], OUT, genes=drugs["genes"],
                                         targets=drugs["targets"], reversal=drugs["reversal"],
                                         cache_dir=CACHE, fdr=STATE["fdr"], reps=reps)
    drugs.pop("proximity_error", None)
    drugs.update(res)
    prox, info = res["proximity"], res["proximity_info"]
    clin = prox[prox["clinical"].astype(bool)] if len(prox) else prox
    return dumps({"info": info, "n_drugs": int(len(prox)), "n_clinical": int(len(clin)),
                  "n_close": int((clin["z"] < -2).sum()) if len(clin) else 0,
                  "proximity": table(prox, "drugs_proximity")})


# ------------------------------------------------------------ step 8: DepMap dependencies

async def run_depmap(lineage: str) -> str:
    drugs = STATE.get("drugs")
    if drugs is None:
        raise ValueError("run the drug leads step first")
    lineage = (lineage or "").strip()
    if lineage:
        libs = [depmap.DEPENDENCY_LIBRARY, *depmap.TISSUE_LIBRARIES]
        await prefetch(*(treatments.ENRICHR_LIBRARY.format(l) for l in libs), label="DepMap via Enrichr")
    res = depmap.dependency_hypotheses(drugs["genes"], lineage or None, OUT, targets=drugs["targets"],
                                       cache_dir=CACHE)
    drugs.update(res)
    info, dm = res["depmap_info"], res["depmap"]
    out = {"info": info}
    if "skipped" not in info and len(dm):
        out.update(n_sig=int((dm["fdr"] < 0.1).sum()), n_up_needed=int(dm["up_and_needed"].sum()),
                   depmap=table(dm, "depmap_dependencies"))
    return dumps(out)


async def lineages() -> str:
    """Tissues of the CCLE cell lines DepMap covers, with line counts."""
    libs = [depmap.DEPENDENCY_LIBRARY, *depmap.TISSUE_LIBRARIES]
    await prefetch(*(treatments.ENRICHR_LIBRARY.format(l) for l in libs), label="DepMap via Enrichr")
    deps, tissues = depmap.load_dependencies(CACHE), depmap.load_tissues(CACHE)
    counts = pd.Series([tissues[c] for c in deps if c in tissues]).value_counts()
    return dumps([[k, int(v)] for k, v in counts.items()])


# ------------------------------------------------------------ step 9: report

def make_report() -> str:
    ds, case, control = STATE["ds"], STATE["case"], STATE["control"]
    sets = STATE["gene_sets"]
    counts = {src: sum(k.startswith(f"{src}:") for k in sets) for src in ("KEGG", "Reactome")}
    notes = list(ds.notes) + [f"run in the browser (Pyodide); GSEA by NumPy prerank, "
                              f"{STATE.get('permutations', 1000)} permutations"]
    drugs = STATE.get("drugs")
    path = report.write_report(OUT, f"{ds.accession}: {case} vs {control}", case, control, STATE["de"],
                               STATE["gsea"], STATE["ora"], ds.samples, source=f"GEO {ds.accession}",
                               notes=notes, fdr=STATE["fdr"], min_abs_lfc=STATE["min_abs_lfc"],
                               gene_set_counts={k: v for k, v in counts.items() if v}, drugs=drugs,
                               drug_error="" if drugs else "the drug leads step was not run",
                               kind=STATE["kind"])
    return Path(path).read_text()
