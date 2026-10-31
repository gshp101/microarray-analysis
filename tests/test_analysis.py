import numpy as np
import pandas as pd

from microarray_pathways import de, enrichment, pathways


def _toy_expression(seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(200)]
    samples = [f"S{i}" for i in range(12)]
    expr = pd.DataFrame(rng.normal(8, 0.3, (200, 12)), index=genes, columns=samples)
    expr.iloc[:20, :6] += 2  # G0-G19 up in the first six samples (case)
    groups = pd.Series(["case"] * 6 + ["control"] * 6, index=samples)
    return expr, groups


def test_differential_expression_finds_shifted_genes():
    expr, groups = _toy_expression()
    res = de.differential_expression(expr, groups, "case", "control")
    assert set(res.index[:20]) == {f"G{i}" for i in range(20)}
    assert (res.loc["G0", "log2fc"] > 1.5) and res.loc["G0", "fdr"] < 0.01
    assert res.attrs["n_case"] == 6 and res.attrs["n_control"] == 6
    assert len(de.significant(res)) == 20


def test_differential_expression_ignores_unlabeled_samples():
    expr, groups = _toy_expression()
    groups.iloc[0] = np.nan
    res = de.differential_expression(expr, groups, "case", "control")
    assert res.attrs["n_case"] == 5


def test_ora_enriches_the_right_set():
    universe = {f"G{i}" for i in range(200)}
    sets = {"hit": {f"G{i}" for i in range(15)}, "miss": {f"G{i}" for i in range(100, 130)}}
    res = enrichment.ora({f"G{i}" for i in range(12)}, universe, sets).set_index("pathway")
    assert res.loc["hit", "hits"] == 12 and res.loc["hit", "fdr"] < 1e-5
    assert res.loc["miss", "hits"] == 0 and res.loc["miss", "pvalue"] == 1.0


def test_gsea_direction():
    expr, groups = _toy_expression()
    res = de.differential_expression(expr, groups, "case", "control")
    sets = {"up_set": {f"G{i}" for i in range(20)}, "random": {f"G{i}" for i in range(100, 130)}}
    g = enrichment.gsea(res, sets, permutations=200, threads=1).set_index("pathway")
    assert g.loc["up_set", "direction"] == "up" and g.loc["up_set", "fdr"] < 0.05


def test_parse_kegg_keeps_metabolic_and_symbols():
    plist = ("path:hsa00010\tGlycolysis / Gluconeogenesis - Homo sapiens (human)\n"
             "hsa04110\tCell cycle - Homo sapiens (human)\n")
    links = "path:hsa00010\thsa:3101\npath:hsa00010\thsa:999\npath:hsa04110\thsa:1017\n"
    genes = ("hsa:3101\tCDS\t10:1..2\tHK3, HKIII; hexokinase 3\n"
             "hsa:999\tCDS\t1:1..2\tsomething-like protein\n"
             "hsa:1017\tCDS\t12:1..2\tCDK2; cyclin dependent kinase 2\n")
    sets = pathways.parse_kegg(plist, links, genes)
    assert sets == {"KEGG: Glycolysis / Gluconeogenesis": {"HK3"}}


def test_parse_reactome_filters_to_metabolism_subtree():
    gmt = ("Glycolysis\tR-HSA-70171\tHK1\tPKM\n"
           "Cell Cycle\tR-HSA-1640170\tCDK2\n"
           "Mouse thing\tR-MMU-1\tHk1\n")
    rel = "R-HSA-1430728\tR-HSA-71387\nR-HSA-71387\tR-HSA-70171\n"
    sets = pathways.parse_reactome(gmt, rel)
    assert sets == {"Reactome: Glycolysis": {"HK1", "PKM"}}


def test_brite_group_and_signaling_scope():
    brite = ("AMetabolism\nB  Carbohydrate metabolism\nC    00010  Glycolysis / Gluconeogenesis\n"
             "AEnvironmental Information Processing\nB  Signal transduction\nC    04010  MAPK signaling pathway\n"
             "C    04370  VEGF signaling pathway\nB  Signaling molecules and interaction\n"
             "C    04080  Neuroactive ligand-receptor interaction\nACellular Processes\nB  Cell growth and death\n"
             "C    04110  Cell cycle\n")
    assert pathways.parse_kegg_brite_group(brite) == {"04010", "04370"}

    plist = "path:hsa00010\tGlycolysis - Homo sapiens (human)\npath:hsa04370\tVEGF signaling pathway\n"
    links = "path:hsa00010\thsa:3101\npath:hsa04370\thsa:3791\n"
    genes = "hsa:3101\tCDS\t10:1..2\tHK3; hexokinase 3\nhsa:3791\tCDS\t4:1..2\tKDR, VEGFR2; kinase insert\n"
    sig = pathways.parse_kegg(plist, links, genes, keep=lambda pid: pid.endswith("04370"))
    assert sig == {"KEGG: VEGF signaling pathway": {"KDR"}}


def test_parse_reactome_takes_several_roots():
    gmt = "Glycolysis\tR-HSA-70171\tHK1\nSignaling by VEGF\tR-HSA-194138\tKDR\nCell Cycle\tR-HSA-1640170\tCDK2\n"
    rel = "R-HSA-1430728\tR-HSA-70171\nR-HSA-162582\tR-HSA-194138\n"
    sets = pathways.parse_reactome(gmt, rel, (pathways.REACTOME_METABOLISM, pathways.REACTOME_SIGNALING))
    assert set(sets) == {"Reactome: Glycolysis", "Reactome: Signaling by VEGF"}


def test_machinery_fraction():
    assert pathways.machinery_fraction({"RPL3", "RPS6", "PSMA1", "HK1"}) == 0.75
    assert pathways.machinery_fraction({"HK1", "PKM"}) == 0.0
