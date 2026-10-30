import pandas as pd

from microarray_pathways import depmap, report


def _genes():
    return pd.DataFrame({"log2fc": [1.2, 0.8, -1.0], "t": [8.0, 5.0, -6.0], "fdr": [1e-5] * 3,
                         "direction": ["up", "up", "down"],
                         "pathways": [["KEGG: Pyrimidine metabolism"], ["KEGG: Purine metabolism"],
                                      ["KEGG: Fatty acid degradation"]]}, index=["DHODH", "GART", "CPT1A"])


def test_parse_tissues_and_line_names():
    t = depmap.parse_tissues(["NCIH1299 LUNG", "22RV1 PROSTATE TenPx24", "BROKEN"])
    assert t == {"NCIH1299": "LUNG", "22RV1": "PROSTATE"}
    assert depmap.norm_line("NCI-H1299") == depmap.norm_line("NCI H1299") == "NCIH1299"


def test_lineage_dependencies_flags_up_and_needed():
    tissues = {f"LUNG{i}": "LUNG" for i in range(10)} | {f"SKIN{i}": "SKIN" for i in range(30)}
    deps = {f"LUNG{i}": {"DHODH", "CPT1A"} if i < 8 else set() for i in range(10)}
    deps |= {f"SKIN{i}": {"GART", "CPT1A"} if i < 15 else {"CPT1A"} for i in range(30)}
    deps["UNKNOWN"] = {"DHODH"}  # no tissue: ignored
    res = depmap.lineage_dependencies(_genes(), deps, tissues, "lung", {"DHODH": ["Leflunomide"]}).set_index("gene")
    assert res.at["DHODH", "frac_lineage"] == 0.8 and res.at["DHODH", "frac_other"] == 0
    assert res.at["DHODH", "pvalue"] < 0.001 and res.at["DHODH", "up_and_needed"]
    assert not res.at["GART", "up_and_needed"]  # other lines need it, lung lines do not
    assert not res.at["CPT1A", "up_and_needed"]  # went down
    assert res.at["DHODH", "drugs"] == "Leflunomide"
    assert res.index[0] == "DHODH"


def test_opposing_drugs_by_gene():
    targets = pd.DataFrame({"drug": ["Leflunomide", "Other"], "opposes_change": [True, False],
                            "targets": ["DHODH (+1.20)", "GART (+0.80)"]})
    assert depmap.opposing_drugs_by_gene(targets) == {"DHODH": ["Leflunomide"]}


def test_depmap_report_section():
    tissues = {f"L{i}": "LUNG" for i in range(10)} | {f"S{i}": "SKIN" for i in range(10)}
    deps = {f"L{i}": {"DHODH"} for i in range(10)} | {f"S{i}": set() for i in range(10)}
    dm = depmap.lineage_dependencies(_genes(), deps, tissues, "LUNG", {"DHODH": ["Leflunomide"]})
    info = {"lineage": "LUNG", "lineage_lines": 10, "lines": 20, "all_lines": 20}
    html = report.depmap_section({"depmap": dm, "depmap_info": info}, "tumor")
    assert "DHODH" in html and "Leflunomide" in html and "lung" in html
    assert "Skipped" in report.depmap_section({"depmap_info": {"skipped": "no tissue"}}, "AD")
    assert report.depmap_section(None, "x") == ""


def test_no_lineage_skips(tmp_path):
    res = depmap.dependency_hypotheses(_genes(), None, tmp_path)
    assert res["depmap"].empty and "skipped" in res["depmap_info"]
