# microarray-pathways

Python tools that take a public gene expression microarray study (GEO) comparing two
groups, such as Alzheimer's vs control or tumor vs normal, and report which metabolic
pathways are affected.

## Status

- Step 1 (load datasets): done; all four curated datasets load with the group sizes
  reported by the original studies.
- Step 2 (differential expression): done. Welch t-test per gene, log2 fold change,
  Benjamini-Hochberg FDR, volcano plot.
- Step 3 (metabolic pathway enrichment): done. KEGG metabolism + Reactome
  "Metabolism" gene sets, scored by over-representation (ORA, up and down genes
  separately) and preranked GSEA.
- Step 4 (HTML report tool): done. One command produces a self-contained HTML report
  from a curated dataset, any GEO series, or your own expression table.

## Usage

    pip install -e .[dev]
    microarray-pathways list                       # curated datasets
    microarray-pathways load alzheimers_brain      # writes data/GSE5281/
    microarray-pathways load GSE12345 --group case='tumou?r' --group control='normal'

    microarray-pathways report alzheimers_brain    # writes results/GSE5281/report.html

    # any GEO series: give two group rules (regex on sample metadata), case group first
    microarray-pathways report GSE15852 --group tumor='tumou?r|cancer' --group normal='normal'

    # your own data: expression table (first column gene symbols, one column per sample)
    # and a sample sheet with "sample" and "group" columns
    microarray-pathways report --expression expr.tsv --samples samples.tsv \
        --case AD --control control --name my_study

`load` writes `expression.tsv.gz` (log2, one row per gene), `samples.tsv` (metadata and
group label) and `notes.txt`, and warns if group sizes differ from the published study.

## Curated datasets

| key | GEO | comparison |
|---|---|---|
| alzheimers_brain | GSE5281 | Alzheimer's vs control, six brain regions |
| lung_adenocarcinoma | GSE10072 | lung adenocarcinoma vs non-tumor |
| breast_cancer | GSE15852 | breast tumor vs adjacent normal |
| colorectal_cancer | GSE32323 | colorectal cancer vs paired normal |

## Analysis outputs

`report` (alias `analyze`) writes to `results/<name>/`:

| file | contents |
|---|---|
| `genes.tsv` | every gene: log2fc, group means, t, p-value, FDR |
| `volcano.png` | volcano plot (FDR < 0.05 and abs(log2FC) >= 1 highlighted) |
| `pathways_gsea.tsv` / `.png` | preranked GSEA on all genes; NES > 0 means up in the case group |
| `pathways_ora.tsv` | over-representation of significant up and down genes |
| `report.html` | self-contained report: plain-language summary, interactive pathway chart, searchable tables, volcano plot, methods and caveats |
| `summary.txt` | sample counts and how many genes and pathways pass |

Gene sets: KEGG human pathways in the Metabolism category (map numbers below 01100;
the global overview maps are excluded) and Reactome pathways under "Metabolism"
(R-HSA-1430728). Reactome sets made up mostly of ribosomal proteins, proteasome
subunits or nucleoporins (e.g. "Selenocysteine synthesis", "Regulation of ODC") are
dropped by default, since they track general protein machinery rather than
metabolic enzymes.

Caveats: the tests treat samples as independent. GSE15852 and GSE32323 are paired
tumor/normal designs and GSE5281 has several brain regions per donor, so a paired or
mixed model would be somewhat more powerful.
