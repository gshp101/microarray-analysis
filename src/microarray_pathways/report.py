"""Single-file HTML report for one analysis.

The report is self-contained (figures embedded, no external scripts), so it can
be emailed or opened offline. It has:
  - a plain-language summary of the metabolic pathways that shift up or down
  - an interactive bar chart of the top pathways (hover for details)
  - searchable, sortable pathway tables (GSEA and over-representation)
  - the volcano plot and the top changed genes
  - the samples used, methods and caveats
"""

from __future__ import annotations

import base64
import datetime as dt
import html
import re
from pathlib import Path

import pandas as pd

from . import __version__


def _esc(x) -> str:
    return html.escape("" if x is None or (isinstance(x, float) and pd.isna(x)) else str(x))


def _fmt_p(x: float) -> str:
    if pd.isna(x):
        return ""
    if x == 0:
        return "&lt;0.001"
    return f"{x:.2g}" if x < 0.01 else f"{x:.3f}"


def _short(name: str) -> str:
    return name.split(": ", 1)[-1]


def _list(names: list[str]) -> str:
    names = [f"<b>{_esc(_short(n))}</b>" for n in names]
    if len(names) <= 1:
        return "".join(names)
    return ", ".join(names[:-1]) + " and " + names[-1]


# ------------------------------------------------------------------ sections

def findings(gsea: pd.DataFrame, case: str, control: str, fdr: float = 0.05, n: int = 5) -> str:
    """Plain-language summary of the strongest pathway shifts."""
    sig = gsea[gsea["fdr"] < fdr]
    if sig.empty:
        return (f"<p>No metabolic pathway shifted at FDR &lt; {fdr} between <b>{_esc(case)}</b> and "
                f"<b>{_esc(control)}</b>. The tables below list the closest candidates.</p>")
    out = [f"<p><b>{len(sig)}</b> of {len(gsea)} metabolic pathways tested shift between "
           f"<b>{_esc(case)}</b> and <b>{_esc(control)}</b> (GSEA, FDR &lt; {fdr}).</p><ul>"]
    for direction, word in (("down", "Lower"), ("up", "Higher")):
        d = sig[sig["direction"] == direction]
        d = d.reindex(d["nes"].abs().sort_values(ascending=False).index)
        if len(d):
            more = f" (plus {len(d) - n} more)" if len(d) > n else ""
            out.append(f'<li><span class="tag {direction}">{word} in {_esc(case)}</span> '
                       f"{_list(d['pathway'].head(n).tolist())}{more}.</li>")
    out.append("</ul>")
    return "".join(out)


def pathway_chart(gsea: pd.DataFrame, n: int = 20) -> str:
    """Horizontal bar chart of NES for the top pathways, as inline SVG with hover tooltips."""
    top = gsea.dropna(subset=["nes"]).nsmallest(n, "fdr").sort_values("nes", ascending=False)
    if top.empty:
        return ""
    row, label_w, plot_w, pad, fdr_w = 26, 330, 420, 22, 90
    h = row * len(top) + 46
    lim = max(3.0, float(top["nes"].abs().max()) * 1.1)
    x0 = label_w + plot_w / 2
    scale = (plot_w / 2) / lim
    parts = [f'<svg viewBox="0 0 {label_w + plot_w + fdr_w} {h}" role="img" '
             f'aria-label="Top metabolic pathways by GSEA normalized enrichment score" class="chart">']
    for v in range(-int(lim), int(lim) + 1):
        x = x0 + v * scale
        parts.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{pad - 6}" y2="{h - 24}" class="grid"/>'
                     f'<text x="{x:.1f}" y="{h - 8}" class="tick" text-anchor="middle">{v}</text>')
    parts.append(f'<line x1="{x0}" x2="{x0}" y1="{pad - 6}" y2="{h - 24}" class="zero"/>'
                 f'<text x="{label_w + plot_w + fdr_w}" y="{pad - 4}" class="tick" text-anchor="end">FDR</text>')
    for i, r in enumerate(top.itertuples()):
        y = pad + i * row
        w = abs(r.nes) * scale
        x = x0 if r.nes > 0 else x0 - w
        lead = str(r.leading_edge).split(";")[:12]
        tip = (f"{r.pathway}\nNES {r.nes:.2f} · FDR {_fmt_p(r.fdr).replace('&lt;', '<')} · {r.size} genes"
               f"\nleading genes: {', '.join(lead)}")
        name = r.pathway if len(r.pathway) <= 52 else r.pathway[:50] + "…"
        q = _fmt_p(r.fdr)
        parts.append(
            f'<g class="bar {r.direction}"><title>{_esc(tip)}</title>'
            f'<rect x="0" y="{y - 2}" width="{label_w + plot_w + fdr_w}" height="{row}" class="hit"/>'
            f'<text x="{label_w - 10}" y="{y + 14}" class="label" text-anchor="end">{_esc(name)}</text>'
            f'<rect x="{x:.1f}" y="{y + 3}" width="{w:.1f}" height="{row - 8}" rx="3"/>'
            f'<text x="{label_w + plot_w + fdr_w}" y="{y + 14}" class="val" text-anchor="end">{q}</text></g>')
    parts.append("</svg>")
    return "".join(parts)


def _table(df: pd.DataFrame, cols: dict[str, str], table_id: str, numeric: set[str] = frozenset(),
           pvals: set[str] = frozenset()) -> str:
    head = "".join(f'<th data-col="{i}" class="{"num" if c in numeric | pvals else ""}">{_esc(t)}</th>'
                   for i, (c, t) in enumerate(cols.items()))
    rows = []
    for r in df.itertuples(index=False):
        cells = []
        for c in cols:
            v = getattr(r, c)
            if c in pvals:
                cells.append(f'<td class="num" data-v="{v}">{_fmt_p(v)}</td>')
            elif c in numeric:
                cells.append(f'<td class="num" data-v="{v}">{v:.2f}</td>' if pd.notna(v) else "<td></td>")
            elif c == "direction":
                cells.append(f'<td><span class="tag {v}">{_esc(v)}</span></td>')
            elif c in ("genes", "leading_edge"):
                g = [x for x in re.split(r"[;,]\s*", str(v)) if x] if pd.notna(v) else []
                more = f" +{len(g) - 15} more" if len(g) > 15 else ""
                cells.append(f'<td class="genes" title="{_esc(", ".join(g))}">{_esc(", ".join(g[:15]))}'
                             f'<span class="more">{more}</span></td>')
            else:
                cells.append(f"<td>{_esc(v)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (f'<div class="tablewrap"><input class="filter" data-table="{table_id}" type="search" '
            f'placeholder="Filter rows…" aria-label="Filter table">'
            f'<div class="scroll"><table id="{table_id}"><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div></div>')


def _yes(v) -> str:
    return "yes" if bool(v) and not pd.isna(v) else ""


def treatment_section(drugs: dict | None, case: str, error: str = "") -> str:
    """Ranked drug candidates from both routes, framed as research leads."""
    warn = ('<div class="card warn"><b>Research leads, not treatment advice.</b> These lists come from '
            'matching gene expression to drug databases. A drug that acts on a changed enzyme, or whose '
            'signature in cultured cells looks like the opposite of the disease, may still do nothing, or '
            'harm, in patients: tissue, dose, toxicity and the direction of cause and effect are all '
            'untested here. Each lead needs literature review and lab work before anything else.</div>')
    if drugs is None:
        why = f" ({_esc(error)})" if error else " (run without <code>--no-drugs</code> to include it)"
        return f'<h2>Treatment hypotheses</h2><p class="muted">This step was skipped{why}.</p>'
    genes, targets, rev = drugs["genes"], drugs["targets"], drugs["reversal"]
    opp = targets[targets["opposes_change"].astype(bool)] if len(targets) else targets
    n_up = int((genes["direction"] == "up").sum()) if len(genes) else 0
    intro = (f"<p><b>{len(genes)}</b> genes drive the shifted pathways ({n_up} up, {len(genes) - n_up} down in "
             f"{_esc(case)}). KEGG lists <b>{len(targets)}</b> drugs acting on them, of which <b>{len(opp)}</b> "
             f"act in the direction that would push the change back.</p>")
    if len(opp):
        top = opp.head(5)["drug"].tolist()
        intro += f"<p>Strongest target-based leads: {_list(top)}.</p>"
    clin = rev[rev["clinical_drug"].astype(bool)] if len(rev) else rev
    if len(clin):
        intro += (f"<p>Clinical drugs whose L1000 signature best reverses the metabolic profile: "
                  f"{_list(clin.head(5)['compound'].tolist())}.</p>")

    def short_paths(v) -> str:
        p = [_short(x) for x in str(v).split("; ") if x]
        return "; ".join(p[:3]) + (f" (+{len(p) - 3} more)" if len(p) > 3 else "")

    def fmt_targets(df):
        return df.assign(opposes_change=df["opposes_change"].map(_yes), marketed=df["marketed"].map(_yes),
                         pathways=df["pathways"].map(short_paths))

    tcols = {"drug": "Drug", "action": "Action", "targets": "Targets (log2FC)", "marketed": "Marketed",
             "l1000_z_metabolic": "L1000 z", "pathways": "Pathways", "efficacy": "KEGG efficacy"}
    target_html = (_table(fmt_targets(opp.head(40)), tcols, "drugtargets", numeric={"l1000_z_metabolic"})
                   if len(opp) else '<p class="card">No KEGG drug acts against the changed pathway genes.</p>')
    rest = targets[~targets["opposes_change"].astype(bool)] if len(targets) else targets
    rest_html = (f'<details><summary>Other drugs acting on these genes ({len(rest)}: action unclear or '
                 f'same direction as the change)</summary>'
                 f'{_table(fmt_targets(rest.head(60)), {**tcols, "opposes_change": "Opposes"}, "drugother", numeric={"l1000_z_metabolic"})}'
                 f'</details>') if len(rest) else ""

    rcols = {"compound": "Compound", "z_metabolic": "Metabolic z", "percentile_metabolic": "Percentile",
             "z_all": "All-gene z", "reversed_genes": "Pathway genes pushed back"}
    rev_fmt = rev.assign(clinical_drug=rev["clinical_drug"].map(_yes)) if len(rev) else rev
    rev_html = (_table(rev_fmt[rev_fmt["clinical_drug"] == "yes"].head(25), rcols, "revclin",
                       numeric={"z_metabolic", "percentile_metabolic", "z_all"})
                + '<details><summary>All compounds, including research tools</summary>'
                + _table(rev_fmt.head(60), {**rcols, "clinical_drug": "Clinical drug"}, "revall",
                         numeric={"z_metabolic", "percentile_metabolic", "z_all"})
                + "</details>") if len(rev) else '<p class="card">No compound signatures could be scored.</p>'

    return f"""
<h2>Treatment hypotheses</h2>
{warn}
<div class="card" style="margin-top:12px">{intro}</div>

<h3>Drugs that act on the changed enzymes</h3>
<p class="muted">Genes in the leading edge of a pathway with FDR &lt; 0.05 that are themselves changed (FDR &lt; 0.05)
the same way, matched to drug targets in KEGG DRUG. Listed when the drug inhibits a gene that went up or activates
one that went down; ranked by marketed status, then by how strongly its targets changed. "L1000 z" is the same
drug's signature score below, when it was profiled (negative = reverses). Full table: <code>drugs_targets.tsv</code>;
genes: <code>drug_target_genes.tsv</code>.</p>
{target_html}
{rest_html}

<h3>Compounds whose expression signature reverses the disease</h3>
<p class="muted">LINCS L1000 consensus signatures: for ~5,400 compounds, the genes each one raises and lowers in
cultured human cell lines. Metabolic z compares the compound's raised genes with its lowered genes on the
{_esc(case)} vs control t statistic, over metabolic pathway genes only; more negative means the compound pushes
those genes back toward control. All-gene z is the same over every gene. Percentile is the share of all compounds
that reverse less strongly. Genes move in co-regulated groups, so treat z as a ranking, not a p-value. Arrows show
how the compound moves each gene. Full table: <code>drugs_reversal.tsv</code>.</p>
{rev_html}
"""


def proximity_section(drugs: dict | None) -> str:
    """Drugs ranked by how close their targets sit to the pathway genes in the interaction network."""
    if drugs is None:
        return ""
    if "proximity" not in drugs:
        err = drugs.get("proximity_error", "")
        return (f'<h3>Drugs whose targets sit close to the pathway genes in the protein network</h3>'
                f'<p class="muted">Network proximity was skipped{" (" + _esc(err) + ")" if err else ""}.</p>')
    prox, info = drugs["proximity"], drugs.get("proximity_info", {})
    clin = prox[prox["clinical"].astype(bool)] if len(prox) else prox
    near = clin[clin["z"] < -2] if len(clin) else clin
    indirect = near[near["n_targets_in_module"] == 0] if len(near) else near
    intro = (f"<p>{info.get('module_genes', 0)} of the {info.get('pathway_genes', 0)} pathway-driving genes are in "
             f"the interaction network. Of <b>{len(clin)}</b> clinical drugs with targets in the network, "
             f"<b>{len(near)}</b> sit closer to them than chance (z &lt; -2), <b>{len(indirect)}</b> of those "
             f"without hitting a pathway gene directly.</p>")
    if len(indirect):
        intro += f"<p>Closest drugs that act through neighbouring proteins: {_list(indirect.drop_duplicates('targets').head(5)['drug'].tolist())}.</p>"

    cols = {"drug": "Drug", "z": "Proximity z", "d_closest": "Distance", "d_random_mean": "Random distance",
            "n_targets_in_module": "Targets in module", "targets": "Targets (module genes first)"}
    numeric = {"z", "d_closest", "d_random_mean", "l1000_z_metabolic"}
    if "opposes_change" in prox:
        cols["opposes_change"] = "Opposes change (KEGG)"
    if "l1000_z_metabolic" in prox:
        cols["l1000_z_metabolic"] = "L1000 z"

    def fmt(df):
        # drugs with the same targets get the same score; show them as one row
        df = df.copy()
        if "opposes_change" in df:
            df["opposes_change"] = df["opposes_change"].map(_yes)
        groups = df.groupby("targets", sort=False)
        df = df[~df["targets"].duplicated()].copy()
        df["drug"] = df["targets"].map(groups["drug"].agg(lambda s: _short_list(", ".join(s), 4)))
        if "opposes_change" in df:
            df["opposes_change"] = df["targets"].map(groups["opposes_change"].agg(lambda s: _yes((s == "yes").any())))
        if "l1000_z_metabolic" in df:
            df["l1000_z_metabolic"] = df["targets"].map(groups["l1000_z_metabolic"].min())
        df["targets"] = df["targets"].map(lambda v: _short_list(v, 8))
        return df

    if len(clin):
        html = _table(fmt(clin.head(120)).head(40), cols, "proxclin", numeric=numeric)
        if len(indirect):
            html += (f'<details><summary>Close without hitting a pathway gene directly ({len(indirect)})</summary>'
                     f'{_table(fmt(indirect.head(150)).head(60), cols, "proxindirect", numeric=numeric)}</details>')
    else:
        html = '<p class="card">No drug targets could be placed in the network.</p>'
    return f"""
<h3>Drugs whose targets sit close to the pathway genes in the protein network</h3>
<div class="card">{intro}</div>
<p class="muted">Network proximity (Guney et al. 2016). The pathway-driving genes above form the disease module. For
each drug, distance is the average, over its targets, of the number of protein interaction steps to the nearest module
gene (0 = the drug hits a module gene). Random distance is the same for {info.get('reps', 1000)} random protein sets
with the same number of genes and the same connectivity, for both the targets and the module; z compares the two, so
drugs on highly connected proteins are not favoured. Being close says the drug could affect these pathways, not in
which direction. KEGG lists some targets only as a whole family (for example every phosphodiesterase for
"PDE inhibitors"); many related targets make the random distances tight and the z large, so read z together with
the targets. Drugs with identical targets share a row (L1000 z: the most reversing of them). Clinical drugs shown;
full table, including the direction cross-references: <code>drugs_proximity.tsv</code>.</p>
{html}
"""


def depmap_section(drugs: dict | None, case: str) -> str:
    """Pathway genes that tissue-matched DepMap cancer cell lines depend on."""
    if drugs is None or "depmap_info" not in drugs:
        return ""
    info, dm = drugs["depmap_info"], drugs.get("depmap", pd.DataFrame())
    head = "<h3>Pathway genes that matching cancer cell lines need to grow (DepMap)</h3>"
    if "skipped" in info:
        return head + f'<p class="muted">Skipped: {_esc(info["skipped"])}. Give <code>--lineage</code> for cancer data.</p>'
    lin = info["lineage"].lower()
    sig = dm[dm["fdr"] < 0.1] if len(dm) else dm
    upn = dm[dm["up_and_needed"].astype(bool)] if len(dm) else dm
    intro = (f"<p>{info['lineage_lines']} {_esc(lin)} cancer cell lines in DepMap (of {info['lines']} with a known "
             f"tissue). <b>{len(sig)}</b> of the {len(dm)} pathway-driving genes are needed by {_esc(lin)} lines "
             f"more often than by other lines (FDR &lt; 0.1); <b>{len(upn)}</b> went up in {_esc(case)} and are "
             f"needed by at least 1 in 10 {_esc(lin)} lines.</p>")
    lead = upn[upn["drugs"].astype(str).str.len() > 0] if len(upn) else upn
    if len(lead):
        intro += ("<p>Up, needed, and with a KEGG drug acting against them: "
                  + ", ".join(f"<b>{_esc(r.gene)}</b> ({_esc(_short_list(r.drugs, 3))})" for r in lead.head(5).itertuples())
                  + ".</p>")
    cols = {"gene": "Gene", "direction": "Direction", "log2fc": "log2FC", "frac_lineage": f"Share of {lin} lines",
            "frac_other": "Share of other lines", "fdr": "FDR", "up_and_needed": "Up and needed",
            "drugs": "KEGG drugs against the change", "pathways": "Pathways"}

    def fmt(df):
        return df.assign(up_and_needed=df["up_and_needed"].map(_yes),
                         drugs=df["drugs"].fillna("").map(lambda v: _short_list(v, 4)),
                         pathways=df["pathways"].fillna("").map(lambda v: "; ".join(_short(x) for x in str(v).split("; ")[:2])))
    shown = dm[(dm["fdr"] < 0.25) | dm["up_and_needed"].astype(bool)] if len(dm) else dm
    table = (_table(fmt(shown.head(60)), cols, "depmap", numeric={"log2fc", "frac_lineage", "frac_other"},
                    pvals={"fdr"}) if len(shown) else
             '<p class="card">No pathway gene stands out as a dependency of these cell lines.</p>')
    return f"""
{head}
<div class="card">{intro}</div>
<p class="muted">DepMap knocks out each gene with CRISPR in about 1,100 cancer cell lines and records which ones the
cells need to keep growing. The list used here (DepMap 2023 CRISPR gene dependency, via Enrichr) leaves out most
genes every cell needs, so it reflects selective dependencies. Shares compare {_esc(lin)} lines with all other lines
of known tissue (one-sided Fisher test, Benjamini-Hochberg FDR across the pathway genes). A gene that went up and
that these cells depend on is a stronger target lead than one that only changed, but cell lines in a dish differ from
tumors, and a dependency in cells does not mean a drug will work in patients. Shown: FDR &lt; 0.25 or up and
needed. Full table: <code>depmap_dependencies.tsv</code>.</p>
{table}
"""


def _short_list(v, n: int) -> str:
    items = [x for x in str(v).split(", ") if x]
    return ", ".join(items[:n]) + (f" (+{len(items) - n} more)" if len(items) > n else "")


def _img(path: Path, alt: str) -> str:
    if not path.exists():
        return ""
    data = base64.b64encode(path.read_bytes()).decode()
    return f'<img src="data:image/png;base64,{data}" alt="{_esc(alt)}">'


def _samples_table(samples: pd.DataFrame, case: str, control: str) -> str:
    s = samples[samples["group"].isin([case, control])]
    cols = [c for c in ("title", "source_name_ch1") if c in s.columns]
    df = s[cols + ["group"]].copy() if cols else s[["group"]].copy()
    df.insert(0, "sample", df.index)
    df.index.name = None
    return _table(df.sort_values(["group", "sample"]),
                  {"sample": "Sample", **{c: c.replace("_ch1", "").replace("_", " ") for c in cols},
                   "group": "Group"}, "samples")


# ------------------------------------------------------------------- page

CSS = """
:root{--bg:#fcfcfb;--panel:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--line:#e6e5e0;--up:#e34948;--down:#2a78d6;
--upbg:#fbe4e3;--downbg:#e1ecfa;--accent:#2a78d6;color-scheme:light}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#1a1a19;--panel:#222220;--ink:#f0efec;
--ink2:#c3c2b7;--line:#383835;--up:#e66767;--down:#3987e5;--upbg:#4a2322;--downbg:#1d3250;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#1a1a19;--panel:#222220;--ink:#f0efec;--ink2:#c3c2b7;--line:#383835;--up:#e66767;
--down:#3987e5;--upbg:#4a2322;--downbg:#1d3250;color-scheme:dark}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1040px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:26px;margin:0 0 4px;letter-spacing:-.01em}h2{font-size:19px;margin:40px 0 10px}
.sub{color:var(--ink2);margin:0 0 20px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:18px 0}
.stat{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.stat b{display:block;font-size:24px;font-variant-numeric:tabular-nums}.stat span{color:var(--ink2);font-size:13px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:16px 18px}
.card ul{margin:8px 0 0;padding-left:18px}.card li{margin:6px 0}
.tag{display:inline-block;font-size:12px;font-weight:600;padding:1px 7px;border-radius:999px}
.tag.up{background:var(--upbg);color:var(--up)}.tag.down{background:var(--downbg);color:var(--down)}
.chart{width:100%;height:auto;display:block}.chart text{fill:var(--ink);font-size:12px}
.chart .tick,.chart .val{fill:var(--ink2);font-size:11px}.chart .grid{stroke:var(--line)}
.chart .zero{stroke:var(--ink2)}.chart .hit{fill:transparent}.chart .bar:hover .hit{fill:var(--line);opacity:.5}
.chart .up rect:not(.hit){fill:var(--up)}.chart .down rect:not(.hit){fill:var(--down)}
.legend{color:var(--ink2);font-size:13px;margin:6px 0 0}
.legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin:0 4px 0 10px;vertical-align:-1px}
.tablewrap{margin-top:8px}.scroll{max-height:460px;overflow:auto;border:1px solid var(--line);border-radius:10px}
.filter{width:100%;max-width:320px;padding:7px 10px;margin-bottom:8px;border:1px solid var(--line);
border-radius:8px;background:var(--panel);color:var(--ink);font:inherit}
table{border-collapse:collapse;width:100%;font-size:13px;background:var(--panel)}
th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{position:sticky;top:0;background:var(--panel);cursor:pointer;user-select:none;white-space:nowrap}
th:after{content:" ↕";color:var(--ink2);font-size:11px}.num{text-align:right;font-variant-numeric:tabular-nums}
td.genes{color:var(--ink2);font-size:12px;max-width:380px}
img{max-width:100%;height:auto;border-radius:8px;background:#fff}
.two{display:grid;grid-template-columns:minmax(0,3fr) minmax(0,2fr);gap:20px;align-items:start}
.two .card{padding:8px}.more{color:var(--accent);white-space:nowrap}
@media (max-width:760px){.two{grid-template-columns:1fr}}
details{margin-top:10px}summary{cursor:pointer;color:var(--accent)}
.muted{color:var(--ink2);font-size:13px}h3{font-size:16px;margin:26px 0 6px}
.card.warn{border-color:var(--up);background:var(--upbg)}.card>p:first-child{margin-top:0}
.card>p:last-child{margin-bottom:0}code{font-size:13px}
.chartwrap{overflow-x:auto}.chartwrap .chart{min-width:620px}
@media (max-width:500px){.stat b{font-size:19px}h1{font-size:22px}}
"""

JS = """
document.querySelectorAll('.filter').forEach(inp=>inp.addEventListener('input',()=>{
 const q=inp.value.toLowerCase();document.querySelectorAll('#'+inp.dataset.table+' tbody tr')
 .forEach(tr=>tr.style.display=tr.textContent.toLowerCase().includes(q)?'':'none')}));
document.querySelectorAll('th').forEach(th=>th.addEventListener('click',()=>{
 const tb=th.closest('table').tBodies[0],i=+th.dataset.col,asc=th.dataset.asc!=='1';th.dataset.asc=asc?'1':'0';
 const key=td=>td.dataset.v!==undefined?parseFloat(td.dataset.v):td.textContent.toLowerCase();
 [...tb.rows].sort((a,b)=>{const x=key(a.cells[i]),y=key(b.cells[i]);return (x>y?1:x<y?-1:0)*(asc?1:-1)})
 .forEach(r=>tb.appendChild(r))}));
"""


def write_report(
    out_dir: str | Path,
    title: str,
    case: str,
    control: str,
    genes: pd.DataFrame,
    gsea: pd.DataFrame,
    ora: pd.DataFrame,
    samples: pd.DataFrame,
    source: str = "",
    notes: list[str] | None = None,
    fdr: float = 0.05,
    min_abs_lfc: float = 1.0,
    gene_set_counts: dict[str, int] | None = None,
    drugs: dict | None = None,
    drug_error: str = "",
) -> Path:
    """Write report.html into out_dir (which already holds the figures and TSVs)."""
    out = Path(out_dir)
    sig = genes[(genes["fdr"] < fdr) & (genes["log2fc"].abs() >= min_abs_lfc)]
    n_up, n_down = int((sig["log2fc"] > 0).sum()), int((sig["log2fc"] < 0).sum())
    n_path = int((gsea["fdr"] < fdr).sum())

    top_genes = pd.concat([sig[sig["log2fc"] > 0].nsmallest(25, "pvalue"),
                           sig[sig["log2fc"] < 0].nsmallest(25, "pvalue")]).sort_values("log2fc", ascending=False)
    top_genes = top_genes.assign(gene=top_genes.index,
                                 direction=["up" if v > 0 else "down" for v in top_genes["log2fc"]])

    gsea_show = gsea[gsea["fdr"] < 0.25]
    ora_show = ora[ora["fdr"] < 0.25]
    mapping = ("Gene symbols as given in the table." if source.startswith("user table")
               else "Probes mapped to gene symbols with the GEO platform annotation.")
    sets_note = ", ".join(f"{v} {k}" for k, v in (gene_set_counts or {}).items())

    body = f"""
<main>
<h1>{_esc(title)}</h1>
<p class="sub">{_esc(case)} vs {_esc(control)}{' · ' + _esc(source) if source else ''} ·
generated {dt.date.today().isoformat()} by microarray-pathways {__version__}</p>

<div class="stats">
 <div class="stat"><b>{genes.attrs.get('n_case', '')} / {genes.attrs.get('n_control', '')}</b>
   <span>samples ({_esc(case)} / {_esc(control)})</span></div>
 <div class="stat"><b>{len(genes):,}</b><span>genes measured</span></div>
 <div class="stat"><b>{n_up:,} ↑ {n_down:,} ↓</b><span>genes changed (FDR &lt; {fdr}, ≥{2 ** min_abs_lfc:g}-fold)</span></div>
 <div class="stat"><b>{n_path} / {len(gsea)}</b><span>metabolic pathways shifted</span></div>
</div>

<h2>Summary</h2>
<div class="card">{findings(gsea, case, control, fdr)}</div>

<h2>Top metabolic pathways</h2>
<div class="card"><div class="chartwrap">{pathway_chart(gsea)}</div>
<p class="legend">Normalized enrichment score from GSEA. <i style="background:var(--up)"></i>higher in {_esc(case)}
<i style="background:var(--down)"></i>lower in {_esc(case)}. Hover a bar for its leading genes.</p></div>

<h2>All pathway results (GSEA)</h2>
<p class="muted">Pathways with FDR &lt; 0.25, the usual GSEA screening threshold; treat FDR between 0.05 and
0.25 as suggestive. Leading genes are those driving the enrichment. Full table: <code>pathways_gsea.tsv</code>.</p>
{_table(gsea_show, {"pathway": "Pathway", "direction": "Direction", "nes": "NES", "fdr": "FDR",
                    "pvalue": "p", "size": "Genes", "leading_edge": "Leading genes"},
        "gsea", numeric={"nes"}, pvals={"fdr", "pvalue"})}

<h2>Over-representation of changed genes</h2>
<p class="muted">Do the significantly up (or down) genes cluster in a pathway more than chance? Hypergeometric
test, pathways with FDR &lt; 0.25. Full table: <code>pathways_ora.tsv</code>.</p>
{_table(ora_show, {"pathway": "Pathway", "direction": "Genes", "hits": "Hits", "size": "Size",
                   "expected": "Expected", "fdr": "FDR", "genes": "Changed genes"},
        "ora", numeric={"expected"}, pvals={"fdr"}) if len(ora_show) else
  '<p class="card">No pathway passes FDR &lt; 0.25 by over-representation. This is common when '
  'changes are spread across many genes; GSEA above is more sensitive to that.</p>'}

{treatment_section(drugs, case, drug_error)}
{proximity_section(drugs)}
{depmap_section(drugs, case)}

<h2>Genes</h2>
<div class="two">
 <div class="card">{_img(out / "volcano.png", "Volcano plot of log2 fold change against p-value")}</div>
 <div>{_table(top_genes, {"gene": "Gene", "direction": "Direction", "log2fc": "log2FC", "fdr": "FDR"},
              "genes", numeric={"log2fc"}, pvals={"fdr"})}
 <p class="muted">Top 25 up and 25 down by p-value. All genes: <code>genes.tsv</code>.</p></div>
</div>

<h2>Methods</h2>
<div class="card">
<ul>
<li><b>Data.</b> {_esc(source) or 'User supplied expression table'}. {mapping} Where several rows measure
one gene, the one with the highest mean is kept. Values are on a log2 scale.</li>
<li><b>Differential expression.</b> Welch's t-test per gene, log2 fold change as the difference of group means,
Benjamini-Hochberg FDR. Changed genes: FDR &lt; {fdr} and at least {2 ** min_abs_lfc:g}-fold.</li>
<li><b>Pathways.</b> {_esc(sets_note) or 'KEGG and Reactome'} metabolic gene sets (KEGG Metabolism category;
Reactome pathways under Metabolism, minus sets made mostly of ribosomal, proteasome or nuclear pore genes),
10 to 500 measured genes each.</li>
<li><b>GSEA.</b> All genes ranked by t statistic; preranked GSEA (gseapy, 1000 permutations). Positive NES
means the pathway's genes sit toward the {_esc(case)}-high end.</li>
<li><b>Over-representation.</b> Hypergeometric test of up and down genes separately against all measured
genes that belong to a pathway.</li>
<li><b>Treatment hypotheses.</b> Drug targets and actions from KEGG DRUG (action read from KEGG's efficacy text,
so a drug with several targets gets one action for all of them). Signature reversal against the LINCS L1000
chemical perturbation consensus signatures (Enrichr library <code>LINCS_L1000_Chem_Pert_Consensus_Sigs</code>),
scored as the difference in mean t between each compound's up and down genes. "Clinical drug" means KEGG DRUG lists
the compound under an official drug name (INN, USAN, JAN), a pharmacopoeia or a trade name. L1000 profiles come
from cancer and other cell lines, so they fit tumor data better than brain or other tissues.</li>
<li><b>Network proximity.</b> Closest distance from each KEGG DRUG drug's targets to the pathway-driving genes in the
human interactome of Menche et al. (2015; about 13,000 proteins and 137,000 interactions after mapping to gene symbols,
largest connected component), compared with 1000 random target sets and modules drawn from the same degree bins
(at least 100 proteins per bin), as in Guney et al. (2016). The interactome is literature-curated, so well studied
proteins have more known interactions.</li>
<li><b>Cell line dependencies.</b> DepMap 2023 CRISPR gene dependency sets per cell line (Enrichr library
<code>DepMap_CRISPR_GeneDependency_CellLines_2023</code>); tissue of each line from its CCLE name. Pathway-driving
genes tested for being needed more often by lines of the matching tissue (one-sided Fisher exact test).</li>
</ul>
<p><b>Caveats.</b> Expression changes suggest, but do not prove, altered pathway activity: enzyme activity also
depends on protein levels, regulation and metabolite supply. Differences in tissue composition between groups
(for example more fat or immune cells in one group) also show up as pathway shifts. Samples are treated as
independent; paired or repeated designs would gain power from a paired model.</p>
</div>

<details><summary>Samples used ({genes.attrs.get('n_case', 0) + genes.attrs.get('n_control', 0)})</summary>
{_samples_table(samples, case, control)}</details>
{'<details><summary>Processing notes</summary><pre class="muted">' + _esc(chr(10).join(notes)) + '</pre></details>' if notes else ''}
</main>"""

    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f"<title>{_esc(title)}</title><style>{CSS}</style></head><body>{body}"
            f"<script>{JS}</script></body></html>")
    path = out / "report.html"
    path.write_text(page, encoding="utf-8")
    return path
