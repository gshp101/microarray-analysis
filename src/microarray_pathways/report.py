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
.muted{color:var(--ink2);font-size:13px}code{font-size:13px}
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
