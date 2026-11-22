# PathRx in the browser

A static web page that runs the PathRx pipeline with [Pyodide](https://pyodide.org)
(Python compiled to WebAssembly). Search GEO for a human microarray series, label the
samples into two groups, then run each step and look at its tables and plots before
moving on:

| Step | What runs | Data source (from the browser) |
|---|---|---|
| 1 Find a dataset | GEO DataSets search | NCBI E-utilities |
| 2 Load the series | `geo.parse_series_matrix`, `geo.parse_platform_table` | GEO FTP over HTTPS |
| 3 Define the groups | `geo.build_dataset` (filter, log2, probes to genes, group regexes) | none |
| 4 Differential expression | `de.run`, volcano plot | none |
| 5 Pathway gene sets | `pathways.reactome_sets`, `pathways.kegg_sets` | Reactome; KEGG via proxy |
| 6 Score pathways | `enrichment.ora_up_down`, `enrichment.gsea(engine="numpy")` | none |
| 7 Drug leads | `treatments.*` (KEGG DRUG targets, LINCS L1000 reversal) | Enrichr; KEGG via proxy |
| 8 Network proximity | `proximity.proximity_hypotheses` | GitHub raw (interactome); KEGG via proxy |
| 9 Cell line dependencies | `depmap.dependency_hypotheses` | Enrichr |
| 10 Report | `report.write_report`, the same HTML the CLI writes | none |

The page calls the `pathrx` package itself; `pathrx_web.py` only downloads files into
the cache that `geo.fetch` reads and turns results into JSON for the page.

## What differs from the command line

- **KEGG needs a proxy.** `rest.kegg.jp` sends no CORS headers, so browsers refuse to
  read it. Without a proxy the page skips the KEGG pathway sets (Reactome only), the
  KEGG DRUG target route, the clinical-drug flag on L1000 compounds, and network
  proximity, and says so at each step. Two ways to get them back:
  - run it locally: `python web/serve.py`, then open http://localhost:8000/. The server
    forwards `/kegg/` to KEGG and the page finds it on its own;
  - or deploy a tiny proxy and paste its URL under "Data sources and KEGG access". For
    example, a Cloudflare Worker:

    ```js
    export default {
      async fetch(req) {
        const url = new URL(req.url);
        const r = await fetch("https://rest.kegg.jp" + url.pathname + url.search);
        return new Response(r.body, {status: r.status, headers: {
          "Content-Type": "text/plain; charset=utf-8", "Access-Control-Allow-Origin": "*"}});
      }
    };
    ```

  KEGG data is free for academic use; commercial use needs a license. Nothing from KEGG
  is bundled with the page.
- **GSEA uses NumPy instead of gseapy**, which has no Pyodide build. `prerank_numpy`
  follows gseapy's method (same gene-label permutations, NES, nominal p and FDR). On
  GSE14520 it gives identical enrichment scores and leading edges, NES within 0.2, and
  81 of 83 pathways significant in both; on GSE10072 the counts match the CLI exactly.
- Each step defaults to the CLI's settings (1000 GSEA permutations, 1000 proximity
  reference sets); faster settings are offered for a quick look.

Downloads are kept in the browser (IndexedDB), so a second run of the same study does not
download again. "Clear downloads" in the settings panel empties it.

## Build and host

    python web/build.py            # writes web/dist/: index.html, worker.js, pathrx.zip, presets.json
    python web/serve.py            # build + serve at http://localhost:8000/ with the KEGG proxy

`web/dist/` is plain static files. `.github/workflows/pages.yml` publishes it to GitHub
Pages on each push to main, once Settings > Pages > Source is set to "GitHub Actions".

Pyodide loads from jsDelivr (`v0.29.5`). Add `?pyodide=<url>/` to the page address to load
it from somewhere else, for example a local copy served with `serve.py --pyodide DIR`.
