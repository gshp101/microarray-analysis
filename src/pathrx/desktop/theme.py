"""Colors and style sheet for the desktop app, shared with the web app's tokens.
Follows the system's light or dark setting."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QGuiApplication, QPalette

LIGHT = dict(bg="#f6f8f8", surface="#ffffff", ink="#14232b", muted="#566a73", line="#d6e0e3",
             accent="#0f7c84", accent_ink="#ffffff", accent_soft="#e0f1f2", up="#d1453b", down="#2f6fd0",
             warn_bg="#fff6e0", warn_line="#e6c36a", err="#b3261e", err_bg="#fdecea", code_bg="#eef3f4",
             case_bg="#fbe9e7", control_bg="#e8f0fc", alt="#f9fbfb")
DARK = dict(bg="#0d1619", surface="#132128", ink="#e4edef", muted="#93a8ae", line="#24363d",
            accent="#4cc3c9", accent_ink="#06282b", accent_soft="#16343a", up="#f0705f", down="#6aa3f5",
            warn_bg="#2e2611", warn_line="#7a6326", err="#ff8a80", err_bg="#3a1a18", code_bg="#0f1c21",
            case_bg="#3a2220", control_bg="#1b2a40", alt="#162730")

_dark: bool | None = None


def is_dark() -> bool:
    global _dark
    if _dark is None:
        try:
            _dark = QGuiApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark
        except AttributeError:  # Qt < 6.5
            _dark = QGuiApplication.palette().color(QPalette.Window).lightness() < 128
    return _dark


def tokens() -> dict[str, str]:
    return DARK if is_dark() else LIGHT


def apply(app) -> None:
    t = tokens()
    app.setStyle("Fusion")
    pal = QPalette()
    for role, key in [(QPalette.Window, "bg"), (QPalette.Base, "surface"), (QPalette.AlternateBase, "alt"),
                      (QPalette.Text, "ink"), (QPalette.WindowText, "ink"), (QPalette.Button, "surface"),
                      (QPalette.ButtonText, "ink"), (QPalette.Highlight, "accent"),
                      (QPalette.HighlightedText, "accent_ink"), (QPalette.ToolTipBase, "surface"),
                      (QPalette.ToolTipText, "ink"), (QPalette.PlaceholderText, "muted"), (QPalette.Link, "accent")]:
        pal.setColor(role, QColor(t[key]))
    pal.setColor(QPalette.Disabled, QPalette.Text, QColor(t["muted"]))
    pal.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(t["muted"]))
    app.setPalette(pal)
    app.setStyleSheet(STYLE.format(**t))


STYLE = """
QWidget {{ font-size: 10.5pt; }}
QMainWindow, QScrollArea, #page {{ background: {bg}; }}
#sidebar {{ background: {surface}; border: none; border-right: 1px solid {line}; outline: 0; padding: 8px 0; }}
#sidebar::item {{ padding: 9px 14px; border-left: 3px solid transparent; color: {ink}; }}
#sidebar::item:selected {{ background: {accent_soft}; border-left: 3px solid {accent}; color: {ink}; }}
#sidebar::item:disabled {{ color: {muted}; }}
#brand {{ font-size: 15pt; font-weight: 700; color: {ink}; }}
#brand-sub {{ color: {muted}; }}
#card {{ background: {surface}; border: 1px solid {line}; border-radius: 10px; }}
#step-n {{ color: {accent}; font-family: monospace; font-weight: 600; }}
#step-title {{ font-size: 17pt; font-weight: 700; }}
#step-state {{ color: {muted}; font-family: monospace; font-weight: 600; letter-spacing: 1px; }}
#step-state[stale="true"] {{ color: {up}; }}
#muted {{ color: {muted}; }}
#why {{ color: {muted}; }}
#warn {{ background: {warn_bg}; border: 1px solid {warn_line}; border-radius: 6px; padding: 8px 10px; }}
#error {{ background: {err_bg}; border: 1px solid {err}; color: {ink}; border-radius: 6px; padding: 8px 10px; }}
#stat {{ background: {bg}; border: 1px solid {line}; border-radius: 8px; }}
#stat-value {{ font-size: 16pt; font-weight: 700; }}
#stat-value-up {{ font-size: 16pt; font-weight: 700; color: {up}; }}
#stat-value-down {{ font-size: 16pt; font-weight: 700; color: {down}; }}
QPushButton {{ padding: 6px 14px; border: 1px solid {line}; border-radius: 6px; background: {surface}; }}
QPushButton:hover {{ border-color: {accent}; }}
QPushButton:disabled {{ color: {muted}; }}
QPushButton#primary {{ background: {accent}; color: {accent_ink}; border-color: {accent}; font-weight: 600; }}
QPushButton#primary:disabled {{ background: {line}; border-color: {line}; color: {muted}; }}
QPushButton#chip {{ padding: 3px 10px; border-radius: 12px; background: {bg}; }}
QPushButton#chip:hover {{ background: {accent_soft}; }}
QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox {{ padding: 4px 6px; border: 1px solid {line}; border-radius: 5px;
    background: {surface}; }}
QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus {{ border-color: {accent}; }}
QLineEdit[invalid="true"] {{ border-color: {err}; }}
QGroupBox {{ border: 1px solid {line}; border-radius: 8px; margin-top: 12px; padding: 10px; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; font-weight: 600; }}
QGroupBox#case {{ border-left: 4px solid {up}; }}
QGroupBox#control {{ border-left: 4px solid {down}; }}
QTableView {{ gridline-color: {line}; border: 1px solid {line}; border-radius: 4px; }}
QHeaderView::section {{ background: {bg}; border: none; border-bottom: 1px solid {line}; padding: 4px 6px;
    font-weight: 600; }}
QTabWidget::pane {{ border: none; }}
QTabBar::tab {{ padding: 6px 14px; border: none; border-bottom: 2px solid transparent; color: {muted}; }}
QTabBar::tab:selected {{ color: {ink}; border-bottom: 2px solid {accent}; }}
QTreeWidget {{ border: 1px solid {line}; border-radius: 4px; }}
QStatusBar {{ background: {surface}; border-top: 1px solid {line}; }}
QProgressBar {{ border: 1px solid {line}; border-radius: 4px; text-align: center; max-height: 14px; }}
QProgressBar::chunk {{ background: {accent}; border-radius: 3px; }}
code {{ background: {code_bg}; }}
"""
