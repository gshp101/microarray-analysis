"""PathRx desktop app: the web app's step-by-step flow, running the pathrx package natively.

    pathrx-desktop            (after pip install "pathrx[desktop]")
    python -m pathrx.desktop

Search GEO, pick two sample groups with a live preview, then run each pipeline step and
look at its tables and plots before moving on, ending with the HTML report. Everything
is written to the results folder like `pathrx report` writes it.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QThreadPool, QTimer, QUrl
from PySide6.QtGui import QAction, QDesktopServices, QIcon, QKeySequence
from PySide6.QtWidgets import (QApplication, QButtonGroup, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog, QFormLayout, QFrame, QGroupBox, QHBoxLayout,
                               QInputDialog, QLabel, QLineEdit, QListWidget, QListWidgetItem, QMainWindow,
                               QMessageBox, QProgressBar, QPushButton, QRadioButton, QScrollArea, QSplitter,
                               QStackedWidget, QTextBrowser, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget)

from .. import __version__, geo
from . import theme
from .engine import Session, StepResult
from .widgets import Chips, Figure, Stats, TableView, Worker, clear, group_color, note, rich, tabs

ICON = Path(__file__).with_name("icons") / "pathrx.png"
APP = "PathRx"
TUTORIAL_URL = "https://gshp101.github.io/microarray-analysis/tutorial.html"


def tutorial_url() -> QUrl:
    """The step-by-step tutorial: the repository's copy when running from a checkout, else the web app's."""
    local = Path(__file__).resolve().parents[3] / "docs" / "tutorial.html"
    return QUrl.fromLocalFile(str(local)) if local.exists() else QUrl(TUTORIAL_URL)


@dataclass(frozen=True)
class Step:
    key: str
    title: str
    why: str
    after: str  # key of the step this one needs
    locked_note: str


STEPS = [
    Step("search", "Find a dataset",
         "Searches GEO DataSets for human series profiled by microarray. Pick one to load, or type an "
         "accession you already know.", "", ""),
    Step("load", "Load the series",
         "Downloads the series matrix (normalized values and sample descriptions) and the platform annotation "
         "that maps each probe to a gene. Large studies are tens of megabytes; downloads are cached.",
         "search", "Pick a dataset first."),
    Step("groups", "Define the groups",
         "Each group is a regular expression matched, ignoring case, against every sample's title, source and "
         "characteristics. Samples that match neither group, or both, are left out. Click a value below to add "
         "it to the selected group.", "load", "Load a series first."),
    Step("de", "Differential expression",
         "For every gene: log2 fold change (case mean minus control mean), Welch's t-test, and "
         "Benjamini-Hochberg false discovery rate.", "groups", "Define the groups first."),
    Step("sets", "Pathway gene sets",
         "KEGG and Reactome pathways for the chosen scope. Sets made mostly of ribosome, proteasome or nuclear "
         "pore genes are dropped because they track housekeeping, not metabolism.", "de",
         "Run differential expression first."),
    Step("path", "Score pathways",
         "Two methods. GSEA ranks every gene by its t statistic and asks whether a pathway's genes bunch up at "
         "the top or bottom; it catches pathways where many genes move a little. ORA asks whether the "
         "significant up (or down) genes fall in a pathway more often than chance.", "sets",
         "Load gene sets first."),
    Step("drugs", "Drug leads",
         "Takes the genes that drive the shifted pathways and finds existing drugs two ways: drugs acting on "
         "those genes in the direction that would push them back (KEGG DRUG), and compounds whose LINCS L1000 "
         "expression signature is the opposite of the disease. These are research leads, not treatment advice.",
         "path", "Score pathways first."),
    Step("prox", "Network proximity",
         "Ranks KEGG drugs by how close their targets sit to the pathway-driving genes in the human protein "
         "interaction network, against random protein sets of the same size and connectivity (Guney et al. "
         "2016). Negative z means closer than chance.", "drugs", "Run drug leads first."),
    Step("depmap", "Cell line dependencies",
         "For cancer data only. Flags pathway-driving genes that DepMap cancer cell lines of the same tissue "
         "need to grow, more often than other cell lines do. Skip it for other diseases.", "drugs",
         "Run drug leads first."),
    Step("report", "Report",
         "The same self-contained HTML report the command line writes, from the steps run so far. It is saved "
         "in the results folder with every table and figure.", "path", "Score pathways first."),
]
INDEX = {s.key: i for i, s in enumerate(STEPS)}


def dependents(key: str) -> list[str]:
    """Steps that need `key`, directly or not."""
    out, todo = [], [key]
    while todo:
        k = todo.pop()
        for s in STEPS:
            if s.after == k and s.key not in out:
                out.append(s.key)
                todo.append(s.key)
    return out


# ------------------------------------------------------------------ one step's page

class StepPage(QWidget):
    def __init__(self, n: int, step: Step):
        super().__init__()
        self.step = step
        self.setObjectName("page")
        card = QFrame()
        card.setObjectName("card")
        self.card_lay = QVBoxLayout(card)
        self.card_lay.setContentsMargins(26, 22, 26, 22)
        self.card_lay.setSpacing(14)

        head = QHBoxLayout()
        num = QLabel(f"{n:02d}")
        num.setObjectName("step-n")
        title = QLabel(step.title)
        title.setObjectName("step-title")
        self.state = QLabel("")
        self.state.setObjectName("step-state")
        head.addWidget(num)
        head.addSpacing(6)
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(self.state)
        self.card_lay.addLayout(head)

        why = QLabel(step.why)
        why.setObjectName("why")
        why.setWordWrap(True)
        self.card_lay.addWidget(why)
        self.locked = note(step.locked_note, "muted")
        self.card_lay.addWidget(self.locked)

        self.controls = QWidget()
        self.controls_lay = QVBoxLayout(self.controls)
        self.controls_lay.setContentsMargins(0, 0, 0, 0)
        self.controls_lay.setSpacing(10)
        self.card_lay.addWidget(self.controls)

        self.output = QWidget()
        self.output_lay = QVBoxLayout(self.output)
        self.output_lay.setContentsMargins(0, 4, 0, 0)
        self.output_lay.setSpacing(12)
        self.card_lay.addWidget(self.output)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(24, 20, 24, 24)
        outer.addWidget(card)
        outer.addStretch(1)
        self.run_buttons: list[QPushButton] = []

    def set_locked(self, locked: bool):
        self.locked.setVisible(locked and bool(self.step.locked_note))
        self.controls.setEnabled(not locked)
        self.controls.setVisible(not locked)
        self.output.setVisible(not locked)

    def set_state(self, text: str, stale: bool = False):
        self.state.setText(text.upper())
        self.state.setProperty("stale", "true" if stale else "false")
        self.state.style().unpolish(self.state)
        self.state.style().polish(self.state)

    def clear_output(self):
        clear(self.output_lay)

    def add(self, w: QWidget):
        self.output_lay.addWidget(w)
        return w

    def show_result(self, r: StepResult, table_name=None, row_color=None):
        """Standard rendering of a StepResult: numbers, notes, text, figures, tables in tabs."""
        self.clear_output()
        if r.stats:
            self.add(Stats(r.stats))
        for n in r.notes:
            self.add(note(n))
        if r.text:
            self.add(rich(r.text))
        for cap, path in r.figures:
            self.add(Figure(path, cap))
        if r.tables:
            parts = [(t, TableView(df, (table_name or self.step.key) + ("" if len(r.tables) == 1 else
                                                                       f"_{re.sub(r'[^a-z0-9]+', '_', t.lower())}"),
                                   row_color)) for t, df in r.tables]
            if len(parts) == 1:
                h = QLabel(parts[0][0])
                h.setStyleSheet("font-weight:700")
                self.add(h)
            self.add(tabs(parts))

    def show_error(self, msg: str, tb: str):
        self.clear_output()
        lab = note(f"<b>This step failed.</b> {msg}", "error")
        lab.setTextFormat(Qt.RichText)
        lab.setToolTip(tb[-3000:])
        self.add(lab)

    def button(self, text: str, primary: bool = True) -> QPushButton:
        """A button that starts work; all of them are disabled while a step runs."""
        b = QPushButton(text)
        if primary:
            b.setObjectName("primary")
        self.run_buttons.append(b)
        return b


def row(*widgets, stretch=True) -> QWidget:
    w = QWidget()
    lay = QHBoxLayout(w)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(10)
    for x in widgets:
        if isinstance(x, str):
            lab = QLabel(x)
            lab.setObjectName("muted")
            lay.addWidget(lab)
        else:
            lay.addWidget(x)
    if stretch:
        lay.addStretch(1)
    return w


# ------------------------------------------------------------------ settings

class SettingsDialog(QDialog):
    def __init__(self, parent, cache: str, out: str):
        super().__init__(parent)
        self.setWindowTitle("PathRx settings")
        self.cache, self.out = QLineEdit(cache), QLineEdit(out)
        form = QFormLayout(self)
        for label, edit in (("Results folder", self.out), ("Download cache", self.cache)):
            b = QPushButton("Choose…")
            b.clicked.connect(lambda _=False, e=edit: self._pick(e))
            form.addRow(label, row(edit, b, stretch=False))
        hint = QLabel("Each dataset gets its own subfolder of the results folder, as with `pathrx report`. "
                      "The cache is shared with the command line.")
        hint.setObjectName("muted")
        hint.setWordWrap(True)
        form.addRow(hint)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        form.addRow(bb)
        self.resize(620, 0)

    def _pick(self, edit):
        d = QFileDialog.getExistingDirectory(self, "Choose folder", edit.text())
        if d:
            edit.setText(d)


# ------------------------------------------------------------------ main window

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.settings = QSettings(APP, APP)
        cache = self.settings.value("cache", str(geo.DEFAULT_CACHE))
        out = self.settings.value("results", str(Path.home() / "PathRx results"))
        self.session = Session(Path(cache), Path(out), progress=self._progress_from_worker)
        self.pool = QThreadPool.globalInstance()
        self.pool.setMaxThreadCount(1)
        self.running: str | None = None
        self.status: dict[str, str] = {s.key: "locked" for s in STEPS}
        self.preset = None
        self.series_fields = []
        self.report_path: Path | None = None

        self.setWindowTitle(f"{APP} — microarray to pathways")
        if ICON.exists():
            self.setWindowIcon(QIcon(str(ICON)))
        self.resize(1320, 900)

        # sidebar
        side = QWidget()
        side.setObjectName("page")
        sl = QVBoxLayout(side)
        sl.setContentsMargins(0, 0, 0, 0)
        sl.setSpacing(0)
        brand = QWidget()
        bl = QVBoxLayout(brand)
        bl.setContentsMargins(18, 18, 18, 12)
        b1 = QLabel(APP)
        b1.setObjectName("brand")
        b2 = QLabel("From a GEO microarray to affected pathways")
        b2.setObjectName("brand-sub")
        b2.setWordWrap(True)
        bl.addWidget(b1)
        bl.addWidget(b2)
        tut = QLabel("<a href='#'>Read the tutorial</a>")
        tut.setToolTip("The whole pipeline worked through on a liver cancer study (F1)")
        tut.linkActivated.connect(lambda _: self.open_tutorial())
        bl.addSpacing(6)
        bl.addWidget(tut)
        sl.addWidget(brand)
        self.nav = QListWidget()
        self.nav.setObjectName("sidebar")
        for i, s in enumerate(STEPS, 1):
            self.nav.addItem(QListWidgetItem(f"{i:02d}   {s.title}"))
        self.nav.currentRowChanged.connect(self._go)
        sl.addWidget(self.nav, 1)
        side.setMinimumWidth(250)
        side.setMaximumWidth(320)

        # pages
        self.stack = QStackedWidget()
        self.pages: dict[str, StepPage] = {}
        for i, s in enumerate(STEPS, 1):
            page = StepPage(i, s)
            self.pages[s.key] = page
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.NoFrame)
            scroll.setWidget(page)
            self.stack.addWidget(scroll)

        split = QSplitter()
        split.addWidget(side)
        split.addWidget(self.stack)
        split.setStretchFactor(1, 1)
        split.setChildrenCollapsible(False)
        self.setCentralWidget(split)

        # status bar
        self.msg = QLabel("Ready.")
        self.bar = QProgressBar()
        self.bar.setFixedWidth(260)
        self.bar.hide()
        self.statusBar().addWidget(self.msg, 1)
        self.statusBar().addPermanentWidget(self.bar)
        self.where = QLabel()
        self.where.setObjectName("muted")
        self.statusBar().addPermanentWidget(self.where)
        self._show_where()

        self._menus()
        self._build_search()
        self._build_load()
        self._build_groups()
        self._build_de()
        self._build_sets()
        self._build_path()
        self._build_drugs()
        self._build_prox()
        self._build_depmap()
        self._build_report()
        self.set_status("search", "ready")
        self.nav.setCurrentRow(0)

    # ---------------------------------------------------------------- plumbing

    def _menus(self):
        m = self.menuBar().addMenu("&File")
        a = QAction("Open results folder", self)
        a.triggered.connect(self.open_results)
        m.addAction(a)
        a = QAction("Settings…", self)
        a.setShortcut(QKeySequence.Preferences)
        a.triggered.connect(self.edit_settings)
        m.addAction(a)
        m.addSeparator()
        a = QAction("Quit", self)
        a.setShortcut(QKeySequence.Quit)
        a.triggered.connect(self.close)
        m.addAction(a)
        h = self.menuBar().addMenu("&Help")
        a = QAction("Tutorial", self)
        a.setShortcut(QKeySequence.HelpContents)
        a.triggered.connect(self.open_tutorial)
        h.addAction(a)
        a = QAction("About PathRx", self)
        a.triggered.connect(self.about)
        h.addAction(a)
        a = QAction("Source on GitHub", self)
        a.triggered.connect(lambda: QDesktopServices.openUrl(QUrl("https://github.com/gshp101/microarray-analysis")))
        h.addAction(a)

    def _show_where(self):
        self.where.setText(f"Results: {self.session.out_root}")

    def _go(self, i: int):
        if i >= 0:
            self.stack.setCurrentIndex(i)

    def goto(self, key: str):
        self.nav.setCurrentRow(INDEX[key])

    def set_status(self, key: str, status: str, label: str = ""):
        """locked, ready, running, done, stale or failed."""
        self.status[key] = status
        page = self.pages[key]
        page.set_locked(status == "locked")
        page.set_state({"done": label or "done", "stale": "out of date", "running": "running…",
                        "failed": "failed"}.get(status, ""), stale=status == "stale")
        item = self.nav.item(INDEX[key])
        mark = {"done": "  ✓", "stale": "  •", "running": "  …", "failed": "  !"}.get(status, "")
        item.setText(f"{INDEX[key] + 1:02d}   {STEPS[INDEX[key]].title}{mark}")
        flags = item.flags()
        item.setFlags(flags & ~Qt.ItemIsEnabled if status == "locked" else flags | Qt.ItemIsEnabled)

    def invalidate(self, key: str):
        """Re-running `key` resets every step that depends on it."""
        for k in dependents(key):
            self.set_status(k, "locked")
            self.pages[k].clear_output()

    def mark_stale(self, *keys: str):
        for k in keys:
            if self.status[k] == "done":
                self.set_status(k, "stale")

    def unlock(self, *keys: str):
        for k in keys:
            if self.status[k] == "locked":
                self.set_status(k, "ready")

    def _progress_from_worker(self, label, loaded, total):
        sig = getattr(self, "_signals", None)
        if sig is not None:
            sig.progress.emit(label, loaded, total)

    def _on_progress(self, label, loaded, total):
        self.msg.setText(label + ("…" if loaded < 0 else ""))
        self.bar.show()
        if loaded < 0 or total <= 0:
            self.bar.setRange(0, 0)
            if loaded > 0:
                self.msg.setText(f"{label}: {loaded / 1048576:.1f} MB")
        else:
            self.bar.setRange(0, 1000)
            self.bar.setValue(int(1000 * min(loaded / total, 1)))
            if total > 1 << 20:
                self.msg.setText(f"{label}: {loaded / 1048576:.1f} of {total / 1048576:.1f} MB")

    def run(self, key: str, fn, *args, done=None, label: str = ""):
        """Run fn(*args) in the background for step `key`; done(result) on success."""
        if self.running:
            self.msg.setText(f"Wait for “{STEPS[INDEX[self.running]].title}” to finish.")
            return
        self.running = key
        self._prev_status = self.status[key]
        self.set_status(key, "running")
        for p in self.pages.values():
            for b in p.run_buttons:
                b.setEnabled(False)
        w = Worker(fn, *args)
        self._worker, self._signals = w, w.signals
        w.signals.progress.connect(self._on_progress)
        w.signals.done.connect(lambda r: self._finish(key, True, r, done, label))
        w.signals.failed.connect(lambda m, tb: self._finish(key, False, (m, tb), done, label))
        self.msg.setText(f"{STEPS[INDEX[key]].title}…")
        self.bar.setRange(0, 0)
        self.bar.show()
        self.pool.start(w)

    def _finish(self, key, ok, result, done, label):
        self.running = None
        self._worker = self._signals = None
        self.bar.hide()
        for p in self.pages.values():
            for b in p.run_buttons:
                b.setEnabled(True)
        if ok:
            self.msg.setText(f"{STEPS[INDEX[key]].title}: done.")
            self.set_status(key, "done", label)
            if done:
                try:
                    done(result)
                except Exception as e:  # rendering problem: show it on the page
                    import traceback
                    self.set_status(key, "failed")
                    self.pages[key].show_error(f"{type(e).__name__}: {e}", traceback.format_exc())
        else:
            msg, tb = result
            self.msg.setText(f"{STEPS[INDEX[key]].title} failed: {msg}")
            self.set_status(key, "failed")
            self.pages[key].show_error(msg, tb)
            print(tb, file=sys.stderr)

    # ---------------------------------------------------------------- 1 search

    def _build_search(self):
        p = self.pages["search"]
        self.q = QLineEdit(placeholderText="e.g. alzheimer hippocampus, hepatocellular carcinoma")
        self.q.returnPressed.connect(self.search)
        sb = QPushButton("Search GEO")
        sb.setObjectName("primary")
        sb.clicked.connect(self.search)
        p.controls_lay.addWidget(row(self.q, sb, stretch=False))

        lab = QLabel("Studies PathRx has checked")
        lab.setStyleSheet("font-weight:600")
        p.controls_lay.addWidget(lab)
        presets = Chips([(f"{d.name or d.key}  {d.accession}", d) for d in Session.presets()])
        presets.clicked.connect(lambda d: self.choose_series(d.accession, d))
        p.controls_lay.addWidget(presets)

        self.acc = QLineEdit(placeholderText="GSE14520")
        self.acc.setMaximumWidth(180)
        self.acc.returnPressed.connect(self._load_accession)
        lb = QPushButton("Load")
        lb.clicked.connect(self._load_accession)
        p.controls_lay.addWidget(row("Or load an accession", self.acc, lb))

        self.hits = QTreeWidget()
        self.hits.setHeaderLabels(["Accession", "Title", "Samples", "Platform", "Date"])
        self.hits.setRootIsDecorated(False)
        self.hits.setAlternatingRowColors(True)
        self.hits.setMinimumHeight(260)
        self.hits.setColumnWidth(0, 100)
        self.hits.setColumnWidth(1, 520)
        self.hits.itemSelectionChanged.connect(self._show_hit)
        self.hits.itemDoubleClicked.connect(lambda *_: self._use_hit())
        self.hit_summary = QTextBrowser()
        self.hit_summary.setMaximumHeight(140)
        use = QPushButton("Use this series")
        use.setObjectName("primary")
        use.clicked.connect(self._use_hit)
        self.hits_box = QWidget()
        hl = QVBoxLayout(self.hits_box)
        hl.setContentsMargins(0, 0, 0, 0)
        self.hits_info = QLabel()
        self.hits_info.setObjectName("muted")
        hl.addWidget(self.hits_info)
        hl.addWidget(self.hits)
        hl.addWidget(self.hit_summary)
        hl.addWidget(row(use))
        self.hits_box.hide()
        p.output_lay.addWidget(self.hits_box)

    def search(self):
        q = self.q.text().strip()
        if not q:
            return
        p = self.pages["search"]
        self.hits_info.setText("Searching GEO…")
        self.hits_box.show()

        def done(r):
            count, hits = r
            self.hits.clear()
            if not hits:
                self.hits_info.setText(f"No human microarray series match “{q}”. Try fewer or broader words.")
                return
            self.hits_info.setText(f"{count:,} series match; showing the {len(hits)} most recent. "
                                   "Double-click one to load it.")
            for h in hits:
                it = QTreeWidgetItem([h["accession"], h["title"], str(h["samples"]), h["platforms"], h["date"]])
                it.setData(0, Qt.UserRole, h)
                it.setToolTip(1, h["title"])
                self.hits.addTopLevelItem(it)
            self.hits.setCurrentItem(self.hits.topLevelItem(0))
            self.set_status("search", "ready")

        self.run("search", self.session.search_geo, q, done=done)
        p.set_state("")

    def _show_hit(self):
        it = self.hits.currentItem()
        if it:
            h = it.data(0, Qt.UserRole)
            self.hit_summary.setHtml(f"<b>{h['accession']}</b> {h['title']}<p>{h['summary']}</p>")

    def _use_hit(self):
        it = self.hits.currentItem()
        if it:
            acc = it.data(0, Qt.UserRole)["accession"]
            self.choose_series(acc, next((d for d in Session.presets() if d.accession == acc), None))

    def _load_accession(self):
        a = self.acc.text().strip().upper()
        if not re.fullmatch(r"GSE\d+", a):
            self.acc.setProperty("invalid", "true")
            self.acc.style().polish(self.acc)
            self.msg.setText("Type a series accession such as GSE14520.")
            return
        self.acc.setProperty("invalid", "false")
        self.acc.style().polish(self.acc)
        self.choose_series(a, next((d for d in Session.presets() if d.accession == a), None))

    # ---------------------------------------------------------------- 2 load

    def _build_load(self):
        pass  # no controls: loading starts from step 1

    def choose_series(self, acc: str, preset=None):
        if self.running:
            self.msg.setText(f"Wait for “{STEPS[INDEX[self.running]].title}” to finish.")
            return
        self.preset = preset
        self.set_status("search", "done", acc)
        self.invalidate("search")
        self.unlock("load")
        self.goto("load")
        page = self.pages["load"]
        page.clear_output()

        def listed(files):
            self.set_status("load", "ready")
            f = files[0]
            if len(files) > 1:
                want = preset.platform if preset else None
                f = next((x for x in files if x["platform"] == want), None)
                if f is None:
                    names = [x["platform"] or x["file"] for x in files]
                    name, ok = QInputDialog.getItem(self, "Pick a platform",
                                                    f"{acc} was run on {len(files)} platforms. PathRx analyzes one "
                                                    "at a time; pick one.", names, 0, False)
                    if not ok:
                        page.show_error("No platform picked.", "")
                        self.set_status("load", "failed")
                        return
                    f = files[names.index(name)]
            QTimer.singleShot(0, lambda: self.run("load", self.session.load_series, acc, f["url"],
                                                  done=self._loaded, label=f.get("platform", "")))

        self.run("load", self.session.list_platforms, acc, done=listed)

    def _loaded(self, r: StepResult):
        self.pages["load"].show_result(r, table_name=f"{self.session.state['accession']}_samples")
        self.series_fields = r.data["fields"]
        self.unlock("groups")
        self.setup_groups()

    # ---------------------------------------------------------------- 3 groups

    def _build_groups(self):
        p = self.pages["groups"]
        self.fields_box = QWidget()
        self.fields_lay = QVBoxLayout(self.fields_box)
        self.fields_lay.setContentsMargins(0, 0, 0, 0)
        p.controls_lay.addWidget(self.fields_box)

        self.target = QButtonGroup(self)
        boxes = []
        self.g = {}
        for key, title in (("case", "Case group"), ("control", "Control group")):
            box = QGroupBox(title)
            box.setObjectName(key)
            name = QLineEdit(key)
            rx = QLineEdit(placeholderText="regular expression, e.g. tumou?r|cancer")
            pick = QRadioButton("add clicked values here")
            self.target.addButton(pick)
            if key == "case":
                pick.setChecked(True)
            count = QLabel("")
            count.setObjectName("muted")
            form = QFormLayout(box)
            form.addRow("Name", name)
            form.addRow("Matches", rx)
            form.addRow(pick, count)
            rx.textEdited.connect(lambda _t, b=pick: b.setChecked(True))
            self.g[key] = (name, rx, pick, count)
            boxes.append(box)
        p.controls_lay.addWidget(row(*boxes, stretch=False))
        self.g_filter = QLineEdit(placeholderText="optional regex")
        self.g_exclude = QLineEdit(placeholderText="optional regex")
        p.controls_lay.addWidget(row("Keep only samples matching", self.g_filter, "Drop samples matching",
                                     self.g_exclude, stretch=False))
        self.g_counts = QLabel()
        self.g_counts.setTextFormat(Qt.RichText)
        self.apply_btn = p.button("Build gene table")
        self.apply_btn.clicked.connect(self.apply_groups)
        p.controls_lay.addWidget(row(self.apply_btn, self.g_counts))
        self.group_table_box = QVBoxLayout()
        p.controls_lay.addLayout(self.group_table_box)

        self._preview_timer = QTimer(self, singleShot=True, interval=250)
        self._preview_timer.timeout.connect(self.preview_groups)
        for w in (self.g["case"][0], self.g["case"][1], self.g["control"][0], self.g["control"][1],
                  self.g_filter, self.g_exclude):
            w.textChanged.connect(lambda _t: self._preview_timer.start())

    def setup_groups(self):
        clear(self.fields_lay)
        if self.series_fields:
            for col, values in self.series_fields:
                lab = QLabel(col)
                lab.setObjectName("muted")
                self.fields_lay.addWidget(lab)
                chips = Chips([((v if len(v) <= 70 else v[:67] + "…") + f"  ({n})", v) for v, n in values])
                chips.clicked.connect(self._add_value)
                self.fields_lay.addWidget(chips)
        else:
            self.fields_lay.addWidget(note("No metadata field repeats across samples; write the patterns by hand "
                                           "from the sample table.", "muted"))
        p = self.preset
        widgets = [self.g["case"][0], self.g["case"][1], self.g["control"][0], self.g["control"][1],
                   self.g_filter, self.g_exclude]
        values = [p.case, p.rules[p.case], p.control, p.rules[p.control], p.sample_filter or "",
                  p.exclude or ""] if p else ["case", "", "control", "", "", ""]
        for w, v in zip(widgets, values):
            w.blockSignals(True)
            w.setText(v)
            w.blockSignals(False)
        if p and p.lineage:
            self.set_lineage(p.lineage)
        self.pages["groups"].clear_output()
        self.preview_groups()

    def _add_value(self, v: str):
        key = "case" if self.g["case"][2].isChecked() else "control"
        rx = self.g[key][1]
        e = re.escape(v)
        rx.setText(f"{rx.text()}|{e}" if rx.text() else e)

    def _group_args(self):
        return (self.g["case"][0].text(), self.g["case"][1].text(), self.g["control"][0].text(),
                self.g["control"][1].text(), self.g_filter.text(), self.g_exclude.text())

    def preview_groups(self):
        if "raw_samples" not in self.session.state:
            return
        args = self._group_args()
        try:
            status = self.session.preview_groups(*args)
        except (re.error, ValueError) as e:
            self.g_counts.setText(f"<span style='color:{theme.tokens()['err']}'>{e}</span>")
            return
        case, control = args[0].strip() or "case", args[2].strip() or "control"
        counts = status.value_counts().to_dict()
        t = theme.tokens()
        colors = {case: t["up"], control: t["down"]}
        order = [case, control, "unassigned", "matches both", "filtered out"]
        self.g_counts.setText("&nbsp;&nbsp;".join(
            f"<span style='color:{colors.get(k, t['muted'])};font-weight:600'>{k}: {counts[k]}</span>"
            for k in order if counts.get(k)))
        self.g["case"][3].setText(f"{counts.get(case, 0)} samples")
        self.g["control"][3].setText(f"{counts.get(control, 0)} samples")
        meta = self.session.state["sample_meta"]
        df = meta.rename_axis("sample").reset_index()
        df.insert(0, "group", status.reindex(meta.index).to_numpy())
        clear(self.group_table_box)
        self.group_table_box.addWidget(TableView(df, "samples_grouped",
                                                 row_color=lambda r: group_color(r["group"], case, control)))

    def apply_groups(self):
        self.invalidate("groups")
        self.run("groups", self.session.apply_groups, *self._group_args(), done=self._grouped)

    def _grouped(self, r: StepResult):
        self.pages["groups"].show_result(r, table_name="expression")
        self.unlock("de")
        self.goto("de")

    # ---------------------------------------------------------------- 4 DE

    def _build_de(self):
        p = self.pages["de"]
        self.fdr = QDoubleSpinBox(decimals=3, minimum=0.001, maximum=0.5, singleStep=0.01, value=0.05)
        self.lfc = QDoubleSpinBox(decimals=2, minimum=0, maximum=5, singleStep=0.25, value=1.0)
        b = p.button("Run")
        b.clicked.connect(lambda: (self.invalidate("de"),
                                   self.run("de", self.session.run_de, self.fdr.value(), self.lfc.value(),
                                            done=self._de_done)))
        p.controls_lay.addWidget(row("FDR below", self.fdr, "|log2 FC| at least", self.lfc, b))

    def _de_done(self, r):
        self.pages["de"].show_result(r, table_name="genes")
        self.unlock("sets")

    # ---------------------------------------------------------------- 5 gene sets

    def _build_sets(self):
        p = self.pages["sets"]
        self.scope = QComboBox()
        for label, v in (("Metabolic", "metabolic"), ("Signaling", "signaling"), ("Metabolic and signaling", "all")):
            self.scope.addItem(label, v)
        self.use_kegg, self.use_reactome = QCheckBox("KEGG"), QCheckBox("Reactome")
        self.use_kegg.setChecked(True)
        self.use_reactome.setChecked(True)
        b = p.button("Load gene sets")
        b.clicked.connect(self._run_sets)
        p.controls_lay.addWidget(row("Pathways", self.scope, self.use_kegg, self.use_reactome, b))

    def _run_sets(self):
        sources = tuple(s for s, c in (("kegg", self.use_kegg), ("reactome", self.use_reactome)) if c.isChecked())
        self.invalidate("sets")
        self.run("sets", self.session.load_gene_sets, self.scope.currentData(), sources, done=self._sets_done)

    def _sets_done(self, r):
        self.pages["sets"].show_result(r, table_name="gene_sets")
        self.unlock("path")

    # ---------------------------------------------------------------- 6 pathways

    def _build_path(self):
        p = self.pages["path"]
        self.perms = QComboBox()
        self.perms.addItem("1000 (as the command line)", 1000)
        self.perms.addItem("200 (quick look)", 200)
        b = p.button("Run")
        b.clicked.connect(lambda: (self.invalidate("path"),
                                   self.run("path", self.session.run_pathways, self.perms.currentData(),
                                            done=self._path_done)))
        p.controls_lay.addWidget(row("GSEA permutations", self.perms, b))

    def _path_done(self, r):
        self.pages["path"].show_result(r, table_name="pathways")
        self.unlock("drugs", "report")

    # ---------------------------------------------------------------- 7 drugs

    def _build_drugs(self):
        p = self.pages["drugs"]
        b = p.button("Run")
        b.clicked.connect(lambda: (self.invalidate("drugs"), self.mark_stale("report"),
                                   self.run("drugs", self.session.run_drugs, done=self._drugs_done)))
        p.controls_lay.addWidget(row(b))

    def _drugs_done(self, r):
        self.pages["drugs"].show_result(r, table_name="drugs")
        self.unlock("prox", "depmap")

    # ---------------------------------------------------------------- 8 proximity

    def _build_prox(self):
        p = self.pages["prox"]
        self.reps = QComboBox()
        self.reps.addItem("1000 (as the command line)", 1000)
        self.reps.addItem("200 (quick look)", 200)
        b = p.button("Run")
        b.clicked.connect(lambda: (self.mark_stale("report"),
                                   self.run("prox", self.session.run_proximity, self.reps.currentData(),
                                            done=lambda r: self.pages["prox"].show_result(r, "drugs_proximity"))))
        p.controls_lay.addWidget(row("Random reference sets", self.reps, b))

    # ---------------------------------------------------------------- 9 DepMap

    def _build_depmap(self):
        p = self.pages["depmap"]
        self.lineage = QComboBox()
        self.lineage.addItem("Not cancer, skip", "")
        self.lineage.setMinimumWidth(240)
        lt = p.button("List tissues", primary=False)
        lt.clicked.connect(self._list_tissues)
        b = p.button("Run")
        b.clicked.connect(lambda: (self.mark_stale("report"),
                                   self.run("depmap", self.session.run_depmap, self.lineage.currentData(),
                                            done=lambda r: self.pages["depmap"].show_result(r, "depmap_dependencies"))))
        p.controls_lay.addWidget(row("Tissue of matching cell lines", self.lineage, lt, b))

    def set_lineage(self, v: str):
        i = self.lineage.findData(v)
        if i < 0:
            self.lineage.addItem(v.lower(), v)
            i = self.lineage.count() - 1
        self.lineage.setCurrentIndex(i)

    def _list_tissues(self):
        cur = self.lineage.currentData()
        prev = self.status["depmap"]

        def done(items):
            self.lineage.clear()
            self.lineage.addItem("Not cancer, skip", "")
            for t, n in items:
                self.lineage.addItem(f"{t.lower()} ({n} lines)", t)
            if cur:
                self.set_lineage(cur)
            self.set_status("depmap", prev if prev != "running" else "ready")
            self.msg.setText(f"{len(items)} tissues listed. Pick one and run.")

        self.run("depmap", self.session.lineages, done=done)

    # ---------------------------------------------------------------- 10 report

    def _build_report(self):
        p = self.pages["report"]
        b = p.button("Build report")
        b.clicked.connect(lambda: self.run("report", self.session.make_report, done=self._report_done))
        self.open_browser = QPushButton("Open in browser")
        self.open_browser.clicked.connect(lambda: self.report_path and QDesktopServices.openUrl(
            QUrl.fromLocalFile(str(self.report_path))))
        self.save_report = QPushButton("Save a copy…")
        self.save_report.clicked.connect(self._save_report)
        folder = QPushButton("Open results folder")
        folder.clicked.connect(self.open_results)
        for w in (self.open_browser, self.save_report):
            w.setEnabled(False)
        p.controls_lay.addWidget(row(b, self.open_browser, self.save_report, folder))

    def _report_done(self, path: Path):
        self.report_path = Path(path)
        for w in (self.open_browser, self.save_report):
            w.setEnabled(True)
        page = self.pages["report"]
        page.clear_output()
        page.add(note(f"Saved to {self.report_path}", "muted"))
        try:
            from PySide6.QtWebEngineWidgets import QWebEngineView
        except ImportError:
            page.add(note("The built-in viewer needs Qt WebEngine (part of the PySide6 package). Use Open in "
                          "browser to read the report."))
            return
        view = QWebEngineView()
        view.setMinimumHeight(820)
        view.load(QUrl.fromLocalFile(str(self.report_path)))
        page.add(view)

    def _save_report(self):
        if not self.report_path:
            return
        dest, _ = QFileDialog.getSaveFileName(self, "Save report", f"{self.session.state['accession']}_report.html",
                                              "HTML (*.html)")
        if dest:
            Path(dest).write_bytes(self.report_path.read_bytes())

    # ---------------------------------------------------------------- menu actions

    def open_results(self):
        d = self.session.out_root / self.session.state["accession"] if "accession" in self.session.state \
            else self.session.out_root
        d.mkdir(parents=True, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(d)))

    def edit_settings(self):
        dlg = SettingsDialog(self, str(self.session.cache_dir), str(self.session.out_root))
        if dlg.exec():
            self.session.cache_dir = Path(dlg.cache.text()).expanduser()
            self.session.out_root = Path(dlg.out.text()).expanduser()
            self.settings.setValue("cache", str(self.session.cache_dir))
            self.settings.setValue("results", str(self.session.out_root))
            self._show_where()

    def open_tutorial(self):
        """Show the tutorial in its own window (or the browser without Qt WebEngine)."""
        url = tutorial_url()
        try:
            from PySide6.QtWebEngineWidgets import QWebEngineView
        except ImportError:
            QDesktopServices.openUrl(url)
            return
        if getattr(self, "_tutorial", None) is None:
            self._tutorial = QWebEngineView()
            self._tutorial.setWindowTitle("PathRx tutorial")
            self._tutorial.setWindowIcon(self.windowIcon())
            self._tutorial.resize(1100, 900)
            self._tutorial.load(url)
        self._tutorial.show()
        self._tutorial.raise_()
        self._tutorial.activateWindow()

    def about(self):
        QMessageBox.about(self, "About PathRx",
                          f"<h3>PathRx {__version__}</h3><p>Find metabolic pathways and drug leads affected in public "
                          "gene expression microarray studies.</p><p>MIT licensed. Data downloaded from KEGG is "
                          "free for academic use only; commercial use needs a KEGG license.</p>"
                          "<p><a href='https://github.com/gshp101/microarray-analysis'>"
                          "github.com/gshp101/microarray-analysis</a></p>")


def main(argv: list[str] | None = None) -> int:
    app = QApplication.instance() or QApplication(sys.argv if argv is None else argv)
    app.setApplicationName(APP)
    app.setOrganizationName(APP)
    if ICON.exists():
        app.setWindowIcon(QIcon(str(ICON)))
    theme.apply(app)
    win = MainWindow()
    win.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
