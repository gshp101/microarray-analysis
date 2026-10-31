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


def test_log2_transform_floors_ratio_data_at_smallest_positive_value():
    # GSE48350-style per-gene normalized ratios: centred on 1, mostly below 1
    ratios = pd.DataFrame({"a": [0.01, 0.5, 1.0, 1400.0], "b": [0.0, 0.8, 1.2, 2.0]})
    out = geo.log2_transform(ratios)
    assert out.loc[1, "a"] == -1 and out.loc[0, "b"] == np.log2(0.01)
    intensities = pd.DataFrame({"a": [0.0, 100.0, 1000.0]})
    assert geo.log2_transform(intensities)["a"].tolist() == [0.0, np.log2(100), np.log2(1000)]


def test_new_dataset_rules_match_real_metadata():
    # titles and characteristics as they appear in the GSE48350 series matrix
    brain = pd.DataFrame({
        "title": ["hippocampus_female_85_AD_12", "Hippocampus_male_91yrs_indiv10",
                  "Hippocampus_male_45yrs_indiv12"],
        "characteristics_ch1": ["age (yrs): 85", "age (yrs): 91", "age (yrs): 45"],
    })
    spec = DATASETS["alzheimers_aging_brain"]
    assert geo.assign_groups(brain, spec.rules).tolist() == ["alzheimers", "control", "control"]
    kept = geo.sample_text(brain).str.contains(spec.sample_filter, case=False, regex=True)
    assert kept.tolist() == [True, True, False]  # controls under 60 left out

    liver = pd.DataFrame({"title": ["Liver Tumor Tissue LCS-079A", "Liver Non-Tumor Tissue LCS-079B"]})
    assert geo.assign_groups(liver, DATASETS["liver_cancer"].rules).tolist() == ["tumor", "nontumor"]
