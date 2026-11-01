# microarray-pathways

Python tools that take a public gene expression microarray study (GEO) comparing two
groups, such as Alzheimer's vs control or tumor vs normal, and report which metabolic
pathways are affected.

## Status

- Step 1 (load datasets): done; all six curated datasets load with the group sizes
  reported by the original studies.
- Step 2 (differential expression): done. Welch t-test per gene, log2 fold change,
  Benjamini-Hochberg FDR, volcano plot.
- Step 3 (metabolic pathway enrichment): done. KEGG metabolism + Reactome
  "Metabolism" gene sets, scored by over-representation (ORA, up and down genes
  separately) and preranked GSEA.
- Step 4 (HTML report tool): done. One command produces a self-contained HTML report
  from a curated dataset, any GEO series, or your own expression table.
- Step 5 (treatment hypotheses): done. Ranks existing drugs two ways, drug-target
  mapping of the changed pathway enzymes (KEGG DRUG) and expression signature
  reversal (LINCS L1000), and adds them to the report as research leads.
- Step 6 (network proximity): done. Ranks KEGG DRUG drugs by how close their targets sit
  to the pathway-driving genes in the human protein interaction network, against
  degree-matched random sets (Guney et al. 2016).
- Step 7 (cell line dependencies): done. For cancer data, flags pathway-driving genes that
  DepMap cancer cell lines of the same tissue need to grow, more often than other lines.
- Signaling pathways (optional): `report --scope all` (or `signaling`) also scores KEGG
  "Signal transduction" pathways and Reactome "Signal Transduction", so drugs acting on
  signaling proteins such as kinase inhibitors can reach the drug lists. The default stays
  `metabolic`. On the liver set (GSE14520) the VEGF, MAPK and PI3K-Akt pathways do not shift
  at the mRNA level (GSEA FDR 0.37 to 0.84), so the approved HCC kinase inhibitors still don't
  appear in the target list.
- Step 8 (trial benchmark): `validate` checks the drug rankings against drugs with known
  positive and failed phase III (or randomized phase II) trials. It reports whether positives
  outrank failures (AUC), and whether the disease's own dataset ranks its positive drugs higher
  than the other datasets do. A benchmark for liver cancer (HCC) ships in `benchmarks/hcc_trials.tsv`.

      microarray-pathways validate GSE14520 --controls GSE10072,GSE15852,GSE32323,GSE5281,GSE48350 \
          --results results_drugs --results results_proximity

  `roc` draws one ROC curve per route with the same benchmark as ground truth (positive-trial
  drugs vs failed-trial drugs). Drugs a route did not score count as not identified and rank
  last. Each panel title gives the AUC and the detections and false alarms in the top 10% of
  the ranking (`--threshold 90`, a percentile). On GSE14520 no route beats chance
  (AUC 0.40 to 0.51).

      microarray-pathways roc GSE14520 --results results_drugs --results results_proximity

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

    # cancer data: name the tissue of matching DepMap cell lines (curated datasets set it)
    microarray-pathways report GSE12345 --group tumor='tumou?r' --group normal='normal' --lineage PANCREAS

`load` writes `expression.tsv.gz` (log2, one row per gene), `samples.tsv` (metadata and
group label) and `notes.txt`, and warns if group sizes differ from the published study.

## Curated datasets

| key | GEO | comparison |
|---|---|---|
| alzheimers_brain | GSE5281 | Alzheimer's vs control, six brain regions |
| lung_adenocarcinoma | GSE10072 | lung adenocarcinoma vs non-tumor |
| breast_cancer | GSE15852 | breast tumor vs adjacent normal |
| colorectal_cancer | GSE32323 | colorectal cancer vs paired normal |
| alzheimers_aging_brain | GSE48350 | Alzheimer's vs cognitively intact controls aged 60+, four brain regions |
| liver_cancer | GSE14520 | hepatocellular carcinoma vs paired non-tumor liver (GPL3921 cohort) |

GSE48350 controls span ages 20 to 99 while the Alzheimer's cases are 60 or older, so
controls under 60 are left out (93 of 173 kept). Its GEO values are GeneSpring per-gene
normalized ratios centred on 1, so they are log2-transformed with a floor at the smallest
positive value rather than at 1. GSE14520 also has a small second cohort on GPL571; only
the 445-sample GPL3921 cohort (225 tumor, 220 non-tumor, 213 patients with both) is used.

## Analysis outputs

`report` (alias `analyze`) writes to `results/<name>/`:

| file | contents |
|---|---|
| `genes.tsv` | every gene: log2fc, group means, t, p-value, FDR |
| `volcano.png` | volcano plot (FDR < 0.05 and abs(log2FC) >= 1 highlighted) |
| `pathways_gsea.tsv` / `.png` | preranked GSEA on all genes; NES > 0 means up in the case group |
| `pathways_ora.tsv` | over-representation of significant up and down genes |
| `report.html` | self-contained report: plain-language summary, interactive pathway chart, searchable tables, volcano plot, methods and caveats |
| `drug_target_genes.tsv` | genes driving the shifted pathways (leading edge of FDR < 0.05 pathways, changed the same way) |
| `drugs_targets.tsv` | KEGG DRUG drugs acting on those genes, with action and whether it opposes the change |
| `drugs_reversal.tsv` | LINCS L1000 compounds ranked by how strongly their signature reverses the disease |
| `drugs_proximity.tsv` | KEGG DRUG drugs ranked by network proximity z of their targets to the pathway-driving genes |
| `depmap_dependencies.tsv` | pathway-driving genes: share of tissue-matched vs other DepMap cell lines that depend on them |
| `summary.txt` | sample counts and how many genes and pathways pass |

Add `--no-drugs` to skip the treatment and proximity steps (they download about 22 MB the first time).

## Treatment hypotheses

Two independent ways to go from the shifted pathways to existing drugs. Both give
leads to research, not treatment advice.

- **Drug targets.** Genes in the GSEA leading edge of a pathway with FDR < 0.05 that are
  themselves changed (FDR < 0.05) in the same direction are matched to drug targets in
  KEGG DRUG. A drug "opposes the change" if KEGG describes it as an inhibitor or
  antagonist and its target went up, or as an activator or agonist and its target went
  down. Ranked by opposing, then marketed (has a trade name or ATC code), then the summed
  |t| of the targets it hits. Most enzyme drugs are inhibitors, so diseases where pathways
  mostly go down (such as the Alzheimer's set) get few opposing drugs.
- **Signature reversal.** The idea behind the Connectivity Map: the LINCS L1000 consensus
  signatures (via Enrichr) list, for about 5,400 compounds, the genes each one raises and
  lowers in cultured cells. For each compound, the mean case-vs-control t of the genes it
  raises minus that of the genes it lowers, scaled to a z-score; strongly negative means
  it pushes the disease profile back toward control. Ranked on metabolic pathway genes,
  with the all-gene score alongside. Genes move in co-regulated groups, so the z-score is
  a ranking, not a calibrated p-value. L1000 profiles come mostly from cancer cell lines,
  so they fit tumor data better than brain tissue.

Compounds found by both routes are cross-referenced (`l1000_z_metabolic` in the targets
table, `known_target_drug` in the reversal table). Tumor vs adjacent normal comparisons
also pick up differences in tissue make-up (for example fat in breast tissue), and
reversal scores reflect that too.

- **Network proximity.** The pathway-driving genes form the disease module. For each
  drug, the closest distance is the mean, over its targets, of the shortest path in the
  human interactome (Menche et al. 2015, the network used by Guney et al. 2016; fetched
  from the paper's GitHub repository) to the nearest module gene. It is compared with
  1000 random target sets and modules drawn from the same degree bins and reported as a
  z-score; z < -2 means closer than chance. This catches drugs that act on neighbours
  of the changed enzymes (regulators, transporters, signalling partners) rather than on
  the enzymes themselves, but says nothing about direction, so check `opposes_change`
  and `l1000_z_metabolic` in `drugs_proximity.tsv`. The interactome is
  literature-curated, so well studied proteins carry more known interactions.

- **Cell line dependencies (DepMap).** For cancer data with a `--lineage` (curated: LUNG,
  BREAST, LARGE INTESTINE, LIVER; none for the Alzheimer's sets), each pathway-driving gene is
  checked against the DepMap 2023 CRISPR gene dependency sets (via Enrichr, which leaves
  out most genes every cell needs). It reports the share of tissue-matched lines that
  depend on the gene, the share of other lines, and a one-sided Fisher test. "Up and
  needed" marks genes that went up in the disease and that at least 1 in 10 matching
  lines need, more often than other lines; KEGG drugs acting against the change are listed
  beside them. Cell line tissue comes from CCLE names, which cover about 620 of the 1,093
  lines. This is a target-priority signal, not synthetic lethality in the strict sense
  (that needs paired knockouts). The DepMap portal and figshare are not reachable from this
  environment; with them, the continuous Chronos gene effect scores and Model.csv
  annotations would replace the Enrichr sets.

Data sources reachable from this project's environment: KEGG REST (rest.kegg.jp) and
Enrichr (maayanlab.cloud), plus raw.githubusercontent.com for the interactome. ChEMBL, DGIdb, Open Targets and the clue.io / SigCom LINCS
APIs would add interaction types, potencies and cell-line-specific signatures if their
hosts are allowed later.

Gene sets: KEGG human pathways in the Metabolism category (map numbers below 01100;
the global overview maps are excluded) and Reactome pathways under "Metabolism"
(R-HSA-1430728). Reactome sets made up mostly of ribosomal proteins, proteasome
subunits or nucleoporins (e.g. "Selenocysteine synthesis", "Regulation of ODC") are
dropped by default, since they track general protein machinery rather than
metabolic enzymes.

Caveats: the tests treat samples as independent. GSE15852, GSE32323 and GSE14520 are
paired tumor/normal designs and GSE5281 and GSE48350 have several brain regions per donor, so a paired or
mixed model would be somewhat more powerful.
