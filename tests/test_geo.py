import numpy as np
import pandas as pd

from microarray_pathways import geo
from microarray_pathways.datasets import DATASETS

MATRIX = """!Series_title\t"toy"
!Series_platform_id\t"GPL0"
!Sample_title\t"AD brain 1"\t"control brain 1"\t"AD brain 2"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"
!Sample_characteristics_ch1\t"disease state: Alzheimer's Disease"\t"disease state: normal"\t"disease state: Alzheimer's Disease"
!Sample_characteristics_ch1\t"region: EC"\t"region: EC"\t"region: HIP"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"
"p1"\t100\t200\t150
"p2"\t50\t60\t55
"p3"\t1000\t900\t950
"p4"\t5\t6\t7
!series_matrix_table_end
"""

PLATFORM = """#ID = probe
ID\tGene title\tGene symbol
p1\tthing\tAAA
p2\tthing\tAAA
p3\tother\tBBB
p4\tambiguous\tCCC /// DDD
"""


def test_parse_series_matrix():
    expr, samples, platform = geo.parse_series_matrix(MATRIX)
    assert platform == "GPL0"
    assert list(expr.columns) == ["GSM1", "GSM2", "GSM3"]
    assert expr.loc["p3", "GSM2"] == 900
    assert samples.loc["GSM3", "characteristics_ch1_2"] == "region: HIP"


def test_platform_and_collapse():
    expr, _, _ = geo.parse_series_matrix(MATRIX)
    sym = geo.parse_platform_table(PLATFORM)
    genes = geo.collapse_to_genes(expr, sym)
    assert list(genes.index) == ["AAA", "BBB"]  # ambiguous probe dropped
    assert genes.loc["AAA", "GSM1"] == 100  # p1 has higher mean than p2


def test_gene_assignment_column():
    text = "ID\tgene_assignment\n1\tNM_1 // TP53 // tumor protein /// NM_2 // X\n2\t---\n"
    sym = geo.parse_platform_table(text)
    assert sym["1"] == "TP53" and sym["2"] == ""


def test_needs_log2():
    expr, _, _ = geo.parse_series_matrix(MATRIX)
    assert geo.needs_log2(expr)
    assert not geo.needs_log2(np.log2(expr))


def test_assign_groups():
    _, samples, _ = geo.parse_series_matrix(MATRIX)
    g = geo.assign_groups(samples, DATASETS["alzheimers_brain"].rules)
    assert g.tolist() == ["alzheimers", "control", "alzheimers"]


def test_lung_rules_match_real_gse10072_metadata():
    # titles/sources as they appear in the GSE10072 series matrix
    samples = pd.DataFrame({
        "title": ["Lung Tumor_GT00006", "Normal Lung_GT00006"],
        "source_name_ch1": ["Adenocarcinoma of the Lung", "Normal Lung Tissue"],
    })
    g = geo.assign_groups(samples, DATASETS["lung_adenocarcinoma"].rules)
    assert g.tolist() == ["tumor", "normal"]


def test_prefix_dir():
    assert geo._prefix_dir("GSE5281") == "GSE5nnn"
    assert geo._prefix_dir("GSE123") == "GSEnnn"
    assert geo._prefix_dir("GPL570") == "GPLnnn"
