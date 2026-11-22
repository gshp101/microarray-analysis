"""Assemble the static PathRx web app into web/dist/ (or the folder given).

    python web/build.py [OUT_DIR]

Writes index.html and worker.js as they are, pathrx.zip (the pathrx package plus
the browser glue module pathrx_web.py, which the page unpacks into Pyodide) and
presets.json (the curated datasets from pathrx/datasets.py). The output is plain
static files: host it on GitHub Pages or any web server.
"""

from __future__ import annotations

import json
import shutil
import sys
import zipfile
from pathlib import Path

WEB = Path(__file__).resolve().parent
ROOT = WEB.parent
PKG = ROOT / "src" / "pathrx"


def presets() -> list[dict]:
    # load datasets.py on its own: it needs only the standard library, so the build
    # runs without NumPy or pandas installed (as on the Pages runner)
    import importlib.util
    spec_ = importlib.util.spec_from_file_location("pathrx_datasets", PKG / "datasets.py")
    datasets = importlib.util.module_from_spec(spec_)
    sys.modules[spec_.name] = datasets
    spec_.loader.exec_module(datasets)

    names = {"alzheimers_brain": "Alzheimer's brain", "alzheimers_aging_brain": "Alzheimer's, aging brain",
             "lung_adenocarcinoma": "Lung adenocarcinoma", "breast_cancer": "Breast cancer",
             "colorectal_cancer": "Colorectal cancer", "liver_cancer": "Liver cancer (HCC)"}
    out = []
    for spec in datasets.DATASETS.values():
        out.append({"key": spec.key, "accession": spec.accession,
                    "name": names.get(spec.key, spec.key.replace("_", " ").capitalize()),
                    "description": spec.description, "case": spec.case, "control": spec.control,
                    "rules": spec.rules, "platform": spec.platform, "sample_filter": spec.sample_filter,
                    "exclude": spec.exclude, "lineage": spec.lineage, "expected": spec.expected})
    return out


def build(out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    for name in ("index.html", "worker.js"):
        shutil.copy2(WEB / name, out / name)
    with zipfile.ZipFile(out / "pathrx.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(PKG.rglob("*")):
            if f.is_file() and "__pycache__" not in f.parts and f.suffix in (".py", ".tsv"):
                z.write(f, Path("pathrx") / f.relative_to(PKG))
        z.write(WEB / "pathrx_web.py", "pathrx_web.py")
    (out / "presets.json").write_text(json.dumps(presets(), indent=1))
    (out / ".nojekyll").write_text("")
    return out


if __name__ == "__main__":
    dest = build(Path(sys.argv[1]) if len(sys.argv) > 1 else WEB / "dist")
    print(f"built {dest}")
