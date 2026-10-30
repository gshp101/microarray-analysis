import numpy as np
import pandas as pd

from microarray_pathways import proximity, report


def _grid(n: int = 30) -> proximity.Network:
    """n x n lattice; node 'r_c'. Plus an isolated edge that the LCC step must drop."""
    edges = []
    for r in range(n):
        for c in range(n):
            if r + 1 < n:
                edges.append((f"{r}_{c}", f"{r + 1}_{c}"))
            if c + 1 < n:
                edges.append((f"{r}_{c}", f"{r}_{c + 1}"))
    return proximity.Network(edges + [("X", "Y")])


def test_network_keeps_largest_component_and_measures_distance():
    net = _grid(5)
    assert len(net) == 25 and "X" not in net.index
    src = np.zeros(len(net), dtype=bool)
    src[net.ids(["0_0", "4_4"])] = True
    d = net.distance_to(src)[:, 0]
    assert d[net.index["0_0"]] == 0 and d[net.index["0_2"]] == 2 and d[net.index["2_2"]] == 4
    assert d[net.index["3_4"]] == 1


def test_parse_sif_maps_entrez_ids():
    edges = proximity.parse_sif("1 1 2\n2 1 3\n3 1 99\n", {"1": "A1BG", "2": "A2M", "3": "NAT1"})
    assert edges == [("A1BG", "A2M"), ("A2M", "NAT1")]


def test_degree_bins_hold_whole_degrees_and_min_size():
    degree = np.array([1] * 60 + [2] * 50 + [3] * 120 + [9] * 5)
    bins = proximity.degree_bins(degree, min_bin_size=100)
    assert list(np.unique(bins[:60])) == [0] and list(np.unique(bins[60:110])) == [0]
    assert len(np.unique(bins[110:])) == 1 and bins[110] == 1  # undersized last bin folded in


def test_drug_near_module_scores_below_far_drug():
    net = _grid(30)
    module = [f"{r}_{c}" for r in range(3) for c in range(3)]
    drugs = pd.DataFrame({"drug": ["Near", "Far", "Direct", "Absent"], "kegg_id": ["D1", "D2", "D3", "D4"],
                          "clinical": [True, True, True, False],
                          "targets": [{"3_1", "1_3"}, {"29_29", "28_29"}, {"1_1"}, {"NOT_IN_NET"}]})
    res = proximity.drug_proximity(net, module, drugs, reps=200).set_index("drug")
    assert "Absent" not in res.index
    assert res.at["Direct", "d_closest"] == 0 and res.at["Direct", "n_targets_in_module"] == 1
    assert res.at["Near", "d_closest"] == 1 and res.at["Near", "n_targets_in_module"] == 0
    assert res.at["Near", "z"] < 0 < res.at["Far", "z"] and res.at["Direct", "z"] < 0


def test_proximity_report_section():
    prox = pd.DataFrame({"drug": ["Near", "Direct"], "kegg_id": ["D1", "D3"], "clinical": [True, True],
                         "z": [-3.0, -4.0], "d_closest": [1.0, 0.0], "d_random_mean": [2.5, 2.5],
                         "d_random_sd": [0.5, 0.6], "n_targets": [2, 1], "n_targets_in_module": [0, 1],
                         "targets": ["3_1, 1_3", "1_1"], "opposes_change": [False, True],
                         "l1000_z_metabolic": [np.nan, -2.0]})
    html = report.proximity_section({"proximity": prox, "proximity_info": {"module_genes": 9, "pathway_genes": 10}})
    assert "Near" in html and "without hitting a pathway gene directly" in html
    assert "skipped (ValueError: x)" in report.proximity_section({"proximity_error": "ValueError: x"})
    assert report.proximity_section(None) == ""
