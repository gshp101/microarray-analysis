"""Reusable Qt widgets for the PathRx desktop app: tables of DataFrames, figures,
stat tiles, wrapping rows of chips, and a background worker."""

from __future__ import annotations

import math
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from PySide6.QtCore import (QAbstractTableModel, QModelIndex, QObject, QPoint, QRect, QRunnable,
                            QSize, QSortFilterProxyModel, Qt, Signal)
from PySide6.QtGui import QBrush, QColor, QPixmap
from PySide6.QtWidgets import (QFileDialog, QFrame, QHBoxLayout, QLabel, QLayout, QLineEdit, QPushButton,
                               QSizePolicy, QTableView, QTabWidget, QVBoxLayout, QWidget)

from . import theme


# ------------------------------------------------------------------ background work

class WorkerSignals(QObject):
    done = Signal(object)
    failed = Signal(str, str)  # message, traceback
    progress = Signal(str, int, int)  # label, loaded (-1 = busy), total


class Worker(QRunnable):
    """Runs fn(*args) off the GUI thread; the result comes back through `signals`."""

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn, self.args, self.kwargs = fn, args, kwargs
        self.signals = WorkerSignals()

    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception as e:  # shown on the step that failed
            self.signals.failed.emit(f"{type(e).__name__}: {e}", traceback.format_exc())
        else:
            self.signals.done.emit(result)


# ------------------------------------------------------------------ tables

def fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, (set, frozenset)):
        return ", ".join(sorted(map(str, v)))
    if isinstance(v, (list, tuple)):
        return ", ".join(map(str, v))
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else ""
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}" if abs(int(v)) >= 10000 else str(int(v))
    if isinstance(v, (float, np.floating)):
        if not math.isfinite(v):
            return ""
        if v == 0:
            return "0"
        a = abs(v)
        return f"{v:.2e}" if a < 1e-3 or a >= 1e6 else f"{v:.4g}"
    return str(v)


class DataFrameModel(QAbstractTableModel):
    def __init__(self, df: pd.DataFrame, row_color=None, parent=None):
        super().__init__(parent)
        self.df = df.reset_index(drop=True)
        self.row_color = row_color  # callable(row Series) -> QColor | None
        self._numeric = [pd.api.types.is_numeric_dtype(t) and not pd.api.types.is_bool_dtype(t)
                         for t in self.df.dtypes]

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.df)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else self.df.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        v = self.df.iat[index.row(), index.column()]
        if role == Qt.DisplayRole:
            return fmt(v)
        if role == Qt.ToolTipRole:
            s = fmt(v)
            return s if len(s) > 40 else None
        if role == Qt.UserRole:  # sort key
            if isinstance(v, (float, np.floating)) and not math.isfinite(v):
                return float("inf")
            return v.item() if hasattr(v, "item") else (v if isinstance(v, (int, float)) else str(v))
        if role == Qt.TextAlignmentRole and self._numeric[index.column()]:
            return int(Qt.AlignRight | Qt.AlignVCenter)
        if role == Qt.BackgroundRole and self.row_color is not None:
            c = self.row_color(self.df.iloc[index.row()])
            return QBrush(c) if c is not None else None
        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        return str(self.df.columns[section]) if orientation == Qt.Horizontal else str(section + 1)


class _Proxy(QSortFilterProxyModel):
    def lessThan(self, left, right):
        a, b = left.data(Qt.UserRole), right.data(Qt.UserRole)
        try:
            return a < b
        except TypeError:
            return str(a) < str(b)


class TableView(QWidget):
    """A DataFrame with a filter box, sortable columns and Save as TSV."""

    def __init__(self, df: pd.DataFrame, name: str = "table", row_color=None, parent=None):
        super().__init__(parent)
        self.df, self.name = df, name
        self.model = DataFrameModel(df, row_color, self)
        self.proxy = _Proxy(self)
        self.proxy.setSourceModel(self.model)
        self.proxy.setFilterKeyColumn(-1)
        self.proxy.setFilterCaseSensitivity(Qt.CaseInsensitive)
        self.view = QTableView()
        self.view.setModel(self.proxy)
        self.view.setSortingEnabled(True)
        self.view.sortByColumn(-1, Qt.AscendingOrder)
        self.view.setAlternatingRowColors(row_color is None)
        self.view.setSelectionBehavior(QTableView.SelectRows)
        self.view.verticalHeader().setDefaultSectionSize(24)
        self.view.horizontalHeader().setStretchLastSection(True)
        self.view.setMinimumHeight(320)
        self._size_columns()

        self.filter = QLineEdit(placeholderText="Filter rows…")
        self.filter.setClearButtonEnabled(True)
        self.filter.textChanged.connect(self._filter)
        self.count = QLabel()
        self.count.setObjectName("muted")
        save = QPushButton("Save as TSV…")
        save.clicked.connect(self.save)
        top = QHBoxLayout()
        top.addWidget(self.filter, 1)
        top.addWidget(self.count)
        top.addWidget(save)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addLayout(top)
        lay.addWidget(self.view)
        self._update_count()

    def _size_columns(self):
        fm, hfm = self.view.fontMetrics(), self.view.horizontalHeader().fontMetrics()
        for i, col in enumerate(self.model.df.columns[:60]):
            sample = [fmt(v) for v in self.model.df.iloc[:50, i]]
            w = max([int(hfm.horizontalAdvance(str(col)) * 1.1) + 12] +
                    [fm.horizontalAdvance(s) for s in sample]) + 24
            self.view.setColumnWidth(i, min(max(w, 60), 360))

    def _filter(self, text):
        self.proxy.setFilterFixedString(text)
        self._update_count()

    def _update_count(self):
        n, total = self.proxy.rowCount(), self.model.rowCount()
        self.count.setText(f"{total:,} rows" if n == total else f"{n:,} of {total:,} rows")

    def save(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save table", f"{self.name}.tsv", "TSV (*.tsv);;CSV (*.csv)")
        if path:
            sep = "," if path.lower().endswith(".csv") else "\t"
            self.df.to_csv(path, sep=sep, index=False, float_format="%.6g")


def tabs(parts: list[tuple[str, QWidget]]) -> QWidget:
    if len(parts) == 1:
        return parts[0][1]
    t = QTabWidget()
    for title, w in parts:
        t.addTab(w, title)
    return t


# ------------------------------------------------------------------ figures and numbers

class Figure(QWidget):
    """A PNG scaled to the page width (never above its own size), with a caption."""

    def __init__(self, path: Path, caption: str = "", parent=None):
        super().__init__(parent)
        self.pix = QPixmap(str(path))
        self.label = QLabel()
        self.label.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.label)
        if caption:
            cap = QLabel(caption)
            cap.setObjectName("muted")
            lay.addWidget(cap)
        self._rescale()

    def _rescale(self):
        if self.pix.isNull():
            self.label.setText("(figure missing)")
            return
        ratio = self.devicePixelRatioF()
        w = min(self.width() or 900, int(self.pix.width() / ratio * 1.0), 1000)
        p = self.pix.scaledToWidth(int(w * ratio), Qt.SmoothTransformation)
        p.setDevicePixelRatio(ratio)
        self.label.setPixmap(p)
        self.label.setFixedHeight(int(p.height() / ratio))

    def resizeEvent(self, e):
        super().resizeEvent(e)
        self._rescale()


class Stats(QWidget):
    """A row of headline numbers."""

    def __init__(self, stats: list[tuple[str, str, str]], parent=None):
        super().__init__(parent)
        lay = FlowLayout(self, spacing=10)
        for value, label, tone in stats:
            tile = QFrame()
            tile.setObjectName("stat")
            v = QLabel(str(value))
            v.setObjectName(f"stat-value-{tone}" if tone else "stat-value")
            lab = QLabel(label)
            lab.setObjectName("muted")
            lab.setWordWrap(True)
            lab.setMinimumWidth(110)
            lab.setMaximumWidth(220)
            tl = QVBoxLayout(tile)
            tl.setContentsMargins(14, 10, 14, 10)
            tl.setSpacing(2)
            tl.addWidget(v)
            tl.addWidget(lab)
            lay.addWidget(tile)


def note(text: str, kind: str = "warn") -> QLabel:
    lab = QLabel(text)
    lab.setObjectName(kind)
    lab.setWordWrap(True)
    lab.setTextInteractionFlags(Qt.TextSelectableByMouse)
    return lab


def rich(text: str) -> QLabel:
    lab = QLabel(text)
    lab.setTextFormat(Qt.RichText)
    lab.setWordWrap(True)
    lab.setOpenExternalLinks(True)
    lab.setTextInteractionFlags(Qt.TextBrowserInteraction)
    return lab


# ------------------------------------------------------------------ chips in wrapping rows

class FlowLayout(QLayout):
    """Lays widgets out left to right, wrapping to new rows (the Qt flow layout example)."""

    def __init__(self, parent=None, spacing=6):
        super().__init__(parent)
        self._items = []
        self._spacing = spacing
        self.setContentsMargins(0, 0, 0, 0)

    def addItem(self, item):
        self._items.append(item)

    def count(self):
        return len(self._items)

    def itemAt(self, i):
        return self._items[i] if 0 <= i < len(self._items) else None

    def takeAt(self, i):
        return self._items.pop(i) if 0 <= i < len(self._items) else None

    def expandingDirections(self):
        return Qt.Orientations(0)

    def hasHeightForWidth(self):
        return True

    def heightForWidth(self, width):
        return self._do_layout(QRect(0, 0, width, 0), True)

    def setGeometry(self, rect):
        super().setGeometry(rect)
        self._do_layout(rect, False)

    def sizeHint(self):
        return self.minimumSize()

    def minimumSize(self):
        size = QSize()
        for item in self._items:
            size = size.expandedTo(item.minimumSize())
        m = self.contentsMargins()
        return size + QSize(m.left() + m.right(), m.top() + m.bottom())

    def _do_layout(self, rect, test):
        x, y, line = rect.x(), rect.y(), 0
        for item in self._items:
            hint = item.sizeHint()
            nx = x + hint.width() + self._spacing
            if nx - self._spacing > rect.right() and line > 0:
                x, y = rect.x(), y + line + self._spacing
                nx, line = x + hint.width() + self._spacing, 0
            if not test:
                item.setGeometry(QRect(QPoint(x, y), hint))
            x, line = nx, max(line, hint.height())
        return y + line - rect.y()


class Chips(QWidget):
    """Wrapping row of small buttons; clicked(value) on press."""

    clicked = Signal(object)

    def __init__(self, items: list[tuple[str, object]], parent=None):
        super().__init__(parent)
        lay = FlowLayout(self)
        for text, value in items:
            b = QPushButton(text)
            b.setObjectName("chip")
            b.setCursor(Qt.PointingHandCursor)
            b.clicked.connect(lambda _=False, v=value: self.clicked.emit(v))
            lay.addWidget(b)


def clear(layout: QLayout) -> None:
    while layout.count():
        item = layout.takeAt(0)
        if item.widget():
            item.widget().deleteLater()
        elif item.layout():
            clear(item.layout())


def group_color(name: str, case: str, control: str):
    if name == case:
        return QColor(theme.tokens()["case_bg"])
    if name == control:
        return QColor(theme.tokens()["control_bg"])
    return None
