import os

import numpy as np
import pandas as pd
import pytest

from pathrx import geo
from pathrx.desktop.engine import Session


def _raw_session(tmp_path, n_genes=300):
    rng = np.random.default_rng(1)
    samples = pd.DataFrame({"title": [f"tumor {i}" for i in range(4)] + [f"normal {i}" for i in range(4)] + ["cell line"],
                            "source_name_ch1": ["liver"] * 9}, index=[f"GSM{i}" for i in range(9)])
    expr = pd.DataFrame(rng.normal(8, 0.3, (n_genes, 9)), index=[f"p{i}" for i in range(n_genes)],
                        columns=samples.index)
    expr.iloc[:40, :4] += 2  # genes G0..G39 up in tumor
    symbols = pd.Series([f"G{i}" for i in range(n_genes)], index=expr.index)
    s = Session(cache_dir=tmp_path / "cache", out_dir=tmp_path / "out")
    s.state = dict(accession="GSE1", url="u", raw_expr=expr, raw_samples=samples, platform="GPL1", symbols=symbols,
                   sample_meta=samples)
    return s


def test_group_status():
    samples = pd.DataFrame({"title": ["tumor a", "normal b", "tumor normal", "other", "cell line tumor"]},
                           index=list("abcde"))
    st = geo.group_status(samples, {"t": "tumor", "n": "normal"}, exclude="cell line")
    assert st.tolist() == ["t", "n", "matches both", "unassigned", "filtered out"]
    assert geo.group_status(samples, {"t": "tumor", "n": ""}).tolist()[:2] == ["unassigned", "unassigned"]
    with pytest.raises(Exception):
        geo.group_status(samples, {"t": "(", "n": "x"})


def test_metadata_fields():
    # titles are all different (more than 12 values), so only the characteristics field is offered
    samples = pd.DataFrame({"title": [f"s{i}" for i in range(15)], "characteristics_ch1": ["x"] * 10 + ["y"] * 5,
                            "other": [1] * 15})
    assert geo.metadata_fields(samples) == [("characteristics_ch1", [("x", 10), ("y", 5)])]


def test_get_retries_transient_errors(monkeypatch):
    calls = []

    class Resp:
        def __init__(self, code):
            self.status_code, self.ok, self.content, self.headers = code, code == 200, b"ok", {}

        def raise_for_status(self):
            if not self.ok:
                raise geo.requests.HTTPError(str(self.status_code))

    def fake_get(url, timeout, stream):
        calls.append(url)
        return Resp(403 if len(calls) < 3 else 200)

    monkeypatch.setattr(geo.requests, "get", fake_get)
    monkeypatch.setattr(geo.time, "sleep", lambda s: None)
    assert geo.get("https://x") == b"ok" and len(calls) == 3


def test_engine_steps_offline(tmp_path):
    s = _raw_session(tmp_path)
    st = s.preview_groups("tumor", "tumor", "normal", "normal", "", "cell line")
    assert st.value_counts().to_dict() == {"tumor": 4, "normal": 4, "filtered out": 1}
    r = s.apply_groups("tumor", "tumor", "normal", "normal", "", "cell line")
    assert r.stats[0][:2] == ("4", "tumor") and (tmp_path / "out/GSE1/distribution_log2.png").exists()
    r = s.run_de(0.05, 1.0)
    assert int(r.stats[1][0]) >= 30 and (tmp_path / "out/GSE1/genes.tsv").exists()
    s.state.update(gene_sets={"Reactome:up set": {f"G{i}" for i in range(30)},
                              "Reactome:null set": {f"G{i}" for i in range(100, 130)},
                              "KEGG:mixed": {f"G{i}" for i in range(20, 50)}},
                   scope="metabolic", kind="metabolic", sources=("reactome",))
    r = s.run_pathways(permutations=100)
    gsea = dict(r.tables)["GSEA"]
    assert gsea.iloc[0]["pathway"] == "Reactome:up set" and gsea.iloc[0]["nes"] > 0
    report = s.make_report()
    assert report.exists() and "up set" in report.read_text()


def test_apply_groups_needs_two_samples(tmp_path):
    s = _raw_session(tmp_path, 20)
    with pytest.raises(ValueError, match="at least 2"):
        s.apply_groups("tumor", "tumor 1", "normal", "normal")


def test_window_flow(tmp_path, monkeypatch):
    pytest.importorskip("PySide6.QtWidgets", exc_type=ImportError)
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    monkeypatch.setenv("HOME", str(tmp_path))
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication

    from pathrx.desktop import app as A

    qa = QApplication.instance() or QApplication([])
    QSettings.setPath(QSettings.NativeFormat, QSettings.UserScope, str(tmp_path))
    w = A.MainWindow()
    assert set(A.dependents("de")) == {"sets", "path", "drugs", "report", "prox", "depmap"}
    assert w.status["search"] == "ready" and w.status["load"] == "locked"

    w.session = _raw_session(tmp_path)
    w.series_fields = geo.metadata_fields(w.session.state["raw_samples"])
    w.unlock("groups")
    w.setup_groups()
    w.g["case"][1].setText("tumor")
    w.g["control"][1].setText("normal")
    w.g_exclude.setText("cell line")
    w.preview_groups()
    assert "case: 4" in w.g_counts.text() and "control: 4" in w.g_counts.text()

    w.set_status("de", "done")
    w.set_status("sets", "done")
    w.invalidate("de")
    assert w.status["sets"] == "locked" and w.status["report"] == "locked"
    w.close()
    qa.processEvents()


def test_tutorial_opens_from_checkout():
    pytest.importorskip("PySide6.QtCore", exc_type=ImportError)
    from pathrx.desktop import app as A

    url = A.tutorial_url()
    assert url.isLocalFile() and url.toLocalFile().endswith("docs/tutorial.html")
