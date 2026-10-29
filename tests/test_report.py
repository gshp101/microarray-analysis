import numpy as np
import pandas as pd

from microarray_pathways import geo, report


def _results():
    genes = pd.DataFrame({"log2fc": [2.0, -1.5, 0.1], "mean_case": [9, 6, 8], "mean_control": [7, 7.5, 8],
                          "t": [8.0, -6.0, 0.3], "pvalue": [1e-8, 1e-6, 0.7], "fdr": [1e-6, 1e-4, 0.8]},
                         index=["HK2", "CPT1A", "ACTB"])
    genes.attrs.update(n_case=3, n_control=3)
    gsea = pd.DataFrame({"pathway": ["KEGG: Glycolysis / Gluconeogenesis", "Reactome: Fatty acid metabolism"],
                         "direction": ["up", "down"], "nes": [2.1, -1.9], "es": [0.6, -0.5],
                         "pvalue": [0.0, 0.002], "fdr": [0.0, 0.01], "size": [60, 150],
                         "leading_edge": ["HK2;PKM;LDHA", "CPT1A;ACADM"]})
    ora = pd.DataFrame({"pathway": ["KEGG: Glycolysis / Gluconeogenesis"], "direction": ["up"], "size": [60],
                        "hits": [1], "expected": [0.1], "pvalue": [0.01], "genes": ["HK2"], "fdr": [0.02]})
    samples = pd.DataFrame({"title": ["a1", "a2", "a3", "c1", "c2", "c3"],
                            "group": ["tumor"] * 3 + ["normal"] * 3}, index=pd.Index([f"S{i}" for i in range(6)], name="sample"))
    return genes, gsea, ora, samples


def test_findings_names_directions():
    _, gsea, _, _ = _results()
    text = report.findings(gsea, "tumor", "normal")
    assert "Higher in tumor" in text and "Glycolysis / Gluconeogenesis" in text
    assert "Lower in tumor" in text and "Fatty acid metabolism" in text


def test_findings_when_nothing_significant():
    _, gsea, _, _ = _results()
    assert "No metabolic pathway" in report.findings(gsea.assign(fdr=0.5), "tumor", "normal")


def test_write_report_is_self_contained(tmp_path):
    genes, gsea, ora, samples = _results()
    path = report.write_report(tmp_path, "Toy <study>", "tumor", "normal", genes, gsea, ora, samples,
                               source="GEO GSE0")
    page = path.read_text()
    assert page.startswith("<!doctype html>")
    assert "Toy &lt;study&gt;" in page  # escaped
    assert "<svg" in page and "CPT1A" in page
    assert "src=\"http" not in page  # no external resources


def test_load_table(tmp_path):
    expr = pd.DataFrame({"gene": ["HK2", "HK2", "PKM"], "s1": [100, 300, 50], "s2": [120, 310, 60],
                         "s3": [90, 280, 55]})
    expr.to_csv(tmp_path / "expr.csv", index=False)
    pd.DataFrame({"sample": ["s1", "s2", "s3"], "group": ["a", "a", "b"]}).to_csv(
        tmp_path / "samples.tsv", sep="\t", index=False)
    ds = geo.load_table(tmp_path / "expr.csv", tmp_path / "samples.tsv")
    assert ds.accession == "expr"
    assert list(ds.expression.index) == ["HK2", "PKM"]
    assert np.isclose(ds.expression.loc["HK2", "s1"], np.log2(300))  # highest-mean row kept, log2 applied
    assert ds.samples["group"].tolist() == ["a", "a", "b"]
