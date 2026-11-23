import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

WEB = Path(__file__).resolve().parents[1] / "web"


def test_build_bundles_package_glue_and_presets(tmp_path):
    sys.path.insert(0, str(WEB))
    import build

    out = build.build(tmp_path)
    names = zipfile.ZipFile(out / "pathrx.zip").namelist()
    assert "pathrx_web.py" in names and "pathrx/geo.py" in names and "pathrx/enrichment.py" in names
    presets = json.loads((out / "presets.json").read_text())
    liver = next(p for p in presets if p["accession"] == "GSE14520")
    assert liver["platform"] == "GPL3921" and liver["case"] in liver["rules"]
    assert (out / "index.html").exists() and (out / "worker.js").exists()
    assert (out / "tutorial.html").exists() and 'href="tutorial.html"' in (out / "index.html").read_text()


def test_group_preview_and_build(monkeypatch):
    sys.path.insert(0, str(WEB))
    import pathrx_web as web

    samples = pd.DataFrame({"title": ["tumor 1", "tumor 2", "normal 1", "normal 2", "cell line"],
                            "source_name_ch1": ["liver"] * 5}, index=[f"GSM{i}" for i in range(5)])
    rng = np.random.default_rng(0)
    expr = pd.DataFrame(rng.normal(8, 1, (4, 5)), index=["p1", "p2", "p3", "p4"], columns=samples.index)
    symbols = pd.Series(["A", "B", "B", ""], index=expr.index)
    web.STATE.clear()
    web.STATE.update(accession="GSE1", url="u", raw_expr=expr, raw_samples=samples, platform="GPL1",
                     symbols=symbols)
    r = json.loads(web.preview_groups("tumor", "tumor", "normal", "normal", "", "cell line"))
    assert r["counts"] == {"tumor": 2, "normal": 2, "filtered out": 1}
    r = json.loads(web.apply_groups("tumor", "tumor", "normal", "normal", "", "cell line"))
    assert (r["n_case"], r["n_control"], r["genes"]) == (2, 2, 2)
    assert json.loads(web.preview_groups("a", "(", "b", "x"))["error"]
