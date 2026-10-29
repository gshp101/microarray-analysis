# microarray-pathways

Python tools that take a public gene expression microarray study (GEO) comparing two
groups, such as Alzheimer's vs control or tumor vs normal, and report which metabolic
pathways are affected.

## Status

- Step 1 (load datasets): loader written and unit tested; real downloads need network
  access to NCBI GEO.
- Step 2 (differential expression), step 3 (metabolic pathway enrichment) and
  step 4 (HTML report tool): not started.

## Usage

    pip install -e .[dev]
    microarray-pathways list                       # curated datasets
    microarray-pathways load alzheimers_brain      # writes data/GSE5281/
    microarray-pathways load GSE12345 --group case='tumou?r' --group control='normal'

`load` writes `expression.tsv.gz` (log2, one row per gene), `samples.tsv` (metadata and
group label) and `notes.txt`, and warns if group sizes differ from the published study.

## Curated datasets

| key | GEO | comparison |
|---|---|---|
| alzheimers_brain | GSE5281 | Alzheimer's vs control, six brain regions |
| lung_adenocarcinoma | GSE10072 | lung adenocarcinoma vs non-tumor |
| breast_cancer | GSE15852 | breast tumor vs adjacent normal |
| colorectal_cancer | GSE32323 | colorectal cancer vs paired normal |
