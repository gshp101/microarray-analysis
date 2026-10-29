import numpy as np
import pandas as pd

from microarray_pathways import report, treatments

KEGG_ENTRY = """ENTRY       D06402                      Drug
NAME        Pemetrexed sodium hydrate (JAN);
            Alimta (TN)
EFFICACY    Antineoplastic, Thymidylate synthase inhibitor, Dihydrofolate reductase inhibitor
CLASS       Antineoplastic
             DG01644  Thymidylate synthase inhibitor
REMARK      ATC code: L01BA04
TARGET      TYMS [HSA:7298] [KO:K00560]
            DHFR [HSA:1719] [KO:K00287]
            GART [HSA:2618] [KO:K11787]
///"""


def _de():
    genes = ["TYMS", "DHFR", "GART", "CPT1A", "ACADM", "HMGCR"] + [f"G{i}" for i in range(40)]
    t = [9.0, 6.0, 5.0, -7.0, -5.0, 0.5] + list(np.linspace(-1, 1, 40))
    de = pd.DataFrame({"t": t, "log2fc": np.array(t) / 4, "fdr": [1e-6] * 5 + [0.5] * 41}, index=genes)
    return de


def _gsea():
    return pd.DataFrame({"pathway": ["KEGG: One carbon pool by folate", "KEGG: Fatty acid degradation"],
                         "direction": ["up", "down"], "nes": [2.2, -2.0], "fdr": [0.001, 0.01],
                         "leading_edge": ["TYMS;DHFR;GART;HMGCR", "CPT1A;ACADM"]})


def test_parse_kegg_drug():
    d = treatments.parse_kegg_drug(KEGG_ENTRY)
    assert d["id"] == "D06402"
    assert d["name"] == "Pemetrexed sodium hydrate"
    assert [t["symbol"] for t in d["targets"]] == ["TYMS", "DHFR", "GART"]
    assert d["marketed"] and d["atc"] == "L01BA04"
    assert treatments.drug_action(d) == "inhibits"


def test_pathway_genes_keeps_changed_leading_edge_genes():
    genes = treatments.pathway_genes(_gsea(), _de())
    assert set(genes.index) == {"TYMS", "DHFR", "GART", "CPT1A", "ACADM"}  # HMGCR not significant
    assert genes.at["CPT1A", "direction"] == "down"


def test_target_candidates_prefers_opposing_drugs():
    genes = treatments.pathway_genes(_gsea(), _de())
    drugs = {"D06402": treatments.parse_kegg_drug(KEGG_ENTRY),
             "D99999": {"name": "Some CPT1 inhibitor", "efficacy": "CPT1 inhibitor", "classes": "",
                        "marketed": True, "atc": ""}}
    res = treatments.target_candidates(genes, {"TYMS": {"D06402"}, "DHFR": {"D06402"}, "CPT1A": {"D99999"}}, drugs)
    assert res.iloc[0]["drug"] == "Pemetrexed"  # salt stripped
    assert bool(res.iloc[0]["opposes_change"])
    assert not res.set_index("drug").at["Some CPT1 inhibitor", "opposes_change"]  # CPT1A already down


def test_signature_library_and_connectivity_sign():
    gmt = ("Reverser Up\t\tCPT1A\tACADM\tG1\tG2\tG3\nReverser Down\t\tTYMS\tDHFR\tGART\tG30\tG31\n"
           "Mimic Up\t\tTYMS\tDHFR\tGART\tG30\tG31\nMimic Down\t\tCPT1A\tACADM\tG1\tG2\tG3\n")
    sigs = treatments.parse_signature_library(gmt)
    assert set(sigs) == {"Reverser", "Mimic"} and sigs["Reverser"]["up"] == {"CPT1A", "ACADM", "G1", "G2", "G3"}
    res = treatments.connectivity(_de()["t"], sigs).set_index("compound")
    assert res.at["Reverser", "z"] < 0 < res.at["Mimic", "z"]
    assert res.at["Reverser", "percentile"] == 100


def test_reversal_candidates_and_report_section():
    gmt = ("Pemetrexed Up\t\tCPT1A\tACADM\tG1\tG2\tG3\nPemetrexed Down\t\tTYMS\tDHFR\tGART\tG30\tG31\n")
    de, gsea = _de(), _gsea()
    genes = treatments.pathway_genes(gsea, de)
    rev = treatments.reversal_candidates(de, treatments.parse_signature_library(gmt), set(de.index), genes,
                                         clinical={"pemetrexed"})
    assert rev.iloc[0]["clinical_drug"] and "TYMS↓" in rev.iloc[0]["reversed_genes"]
    targets = treatments.target_candidates(genes, {"TYMS": {"D06402"}},
                                           {"D06402": treatments.parse_kegg_drug(KEGG_ENTRY)})
    targets["l1000_z_metabolic"] = rev.iloc[0]["z_metabolic"]
    html = report.treatment_section({"genes": genes, "targets": targets, "reversal": rev}, "tumor")
    assert "Research leads, not treatment advice" in html and "Pemetrexed" in html
    assert "skipped" in report.treatment_section(None, "tumor", "ConnectionError: offline")
