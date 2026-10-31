import pandas as pd
import pytest

from microarray_pathways import validation


def _bench(tmp_path):
    path = tmp_path / "bench.tsv"
    pd.DataFrame({"drug": ["Good A", "Good B", "Bad C"],
                  "outcome": ["positive", "positive", "negative"]}).to_csv(path, sep="\t", index=False)
    return validation.load_benchmark(path)


def _write(root, acc, z):
    """L1000 table where lower z ranks higher; z maps compound -> score."""
    (root / acc).mkdir(parents=True)
    pd.DataFrame({"compound": list(z), "z_metabolic": list(z.values())}).to_csv(
        root / acc / "drugs_reversal.tsv", sep="\t", index=False)


def test_benchmark_rejects_unknown_outcomes(tmp_path):
    path = tmp_path / "bench.tsv"
    pd.DataFrame({"drug": ["x"], "outcome": ["maybe"]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError):
        validation.load_benchmark(path)


def test_percentiles_keep_best_row_and_follow_table_order_for_kegg():
    l1000 = pd.DataFrame({"compound": ["A", "a", "B", "C"], "z_metabolic": [1.0, -3.0, 0.0, 2.0]})
    pct = validation.route_percentiles(l1000, "l1000_reversal")
    assert pct["a"] == 100 and pct["c"] == 25
    kegg = pd.DataFrame({"drug": ["first", "second"]})
    assert validation.route_percentiles(kegg, "kegg_targets")["first"] == 100


def test_validate_scores_discrimination_and_target_lift(tmp_path):
    bench = _bench(tmp_path)
    filler = {f"other{i}": float(i) for i in range(10)}
    # target dataset: positives on top, negative at the bottom; control dataset: reversed
    _write(tmp_path, "T", {"Good A": -10.0, "Good B": -9.0, "Bad C": 99.0, **filler})
    _write(tmp_path, "N", {"Good A": 99.0, "Good B": 98.0, "Bad C": -10.0, **filler})
    ranks, summary = validation.validate(bench, ["T", "N"], "T", [tmp_path])
    row = summary.set_index("route").loc["l1000_reversal"]
    assert row["auc_pos_vs_neg"] == 1.0 and row["auc_pos_vs_unlabeled"] == 1.0
    assert row["positives_scored"] == "2/2" and row["median_lift_positive"] > 50
    assert set(ranks["dataset"]) == {"T", "N"}

    out = validation.write(ranks, summary, "T", tmp_path / "out")
    wide = pd.read_csv(out / "benchmark_ranks.tsv", sep="\t")
    assert list(wide.columns[3:]) == ["T", "N"]
