"""Download and parse public GEO microarray series into a gene x sample matrix.

The pipeline for one series is:
  1. download the series matrix (normalized expression + sample metadata)
  2. download the platform annotation and map probe IDs to gene symbols
  3. log2-transform if the values look unlogged (same rule GEO2R uses)
  4. collapse probes to one row per gene (the probe with highest mean expression)
  5. label samples into two groups using regex rules on their metadata
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import requests

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
DEFAULT_CACHE = Path.home() / ".cache" / "microarray_pathways"


@dataclass
class ExpressionDataset:
    accession: str
    platform: str
    expression: pd.DataFrame  # genes x samples, log2 scale
    samples: pd.DataFrame  # one row per sample: metadata text + "group"
    notes: list[str] = field(default_factory=list)

    def group_counts(self) -> pd.Series:
        return self.samples["group"].value_counts(dropna=False)

    def save(self, out_dir: str | Path) -> Path:
        out = Path(out_dir) / self.accession
        out.mkdir(parents=True, exist_ok=True)
        self.expression.to_csv(out / "expression.tsv.gz", sep="\t", float_format="%.4f")
        self.samples.to_csv(out / "samples.tsv", sep="\t")
        (out / "notes.txt").write_text("\n".join([f"platform: {self.platform}", *self.notes]) + "\n")
        return out

    @classmethod
    def read(cls, path: str | Path) -> "ExpressionDataset":
        """Read a dataset written by save()."""
        path = Path(path)
        expression = pd.read_csv(path / "expression.tsv.gz", sep="\t", index_col=0)
        samples = pd.read_csv(path / "samples.tsv", sep="\t", index_col=0, dtype=str)
        lines = (path / "notes.txt").read_text().splitlines()
        platform = lines[0].removeprefix("platform: ") if lines else ""
        return cls(path.name, platform, expression, samples, lines[1:])


# ---------------------------------------------------------------- downloading

def _prefix_dir(acc: str) -> str:
    """GSE5281 -> GSE5nnn, GSE123 -> GSEnnn, GPL570 -> GPLnnn."""
    m = re.fullmatch(r"(GSE|GPL)(\d+)", acc)
    if not m:
        raise ValueError(f"not a GSE/GPL accession: {acc}")
    kind, num = m.groups()
    return f"{kind}{num[:-3]}nnn"


def fetch(url: str, cache_dir: Path = DEFAULT_CACHE, timeout: int = 300) -> bytes:
    """GET a URL, caching the body on disk. Transparently gunzips .gz files."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / re.sub(r"[^A-Za-z0-9._-]", "_", url.split("://", 1)[-1])
    if cached.exists():
        data = cached.read_bytes()
    else:
        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        data = resp.content
        cached.write_bytes(data)
    if url.endswith(".gz") or data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data


def series_matrix_urls(gse: str, cache_dir: Path = DEFAULT_CACHE) -> list[str]:
    """List the series matrix files for a GSE (one per platform)."""
    base = f"{GEO_FTP}/series/{_prefix_dir(gse)}/{gse}/matrix/"
    listing = requests.get(base, timeout=60)
    listing.raise_for_status()
    names = sorted(set(re.findall(r'href="([^"]+_series_matrix\.txt\.gz)"', listing.text)))
    if not names:
        raise FileNotFoundError(f"no series matrix files under {base}")
    return [base + n for n in names]


# ------------------------------------------------------------------- parsing

def parse_series_matrix(text: str) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Parse a GEO series matrix file.

    Returns (expression probes x samples, sample metadata, platform id).
    Repeated metadata keys such as !Sample_characteristics_ch1 are numbered
    characteristics_ch1_1, characteristics_ch1_2, ...
    """
    meta: dict[str, list[str]] = {}
    platform = ""
    table_lines: list[str] = []
    in_table = False
    for line in text.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            table_lines.append(line)
        elif line.startswith("!Sample_"):
            key, *vals = line.split("\t")
            key = key[len("!Sample_"):]
            vals = [v.strip('"') for v in vals]
            n = 1
            name = key
            while name in meta:
                n += 1
                name = f"{key}_{n}"
            meta[name] = vals
        elif line.startswith("!Series_platform_id"):
            platform = line.split("\t", 1)[1].strip().strip('"')

    if not table_lines:
        raise ValueError("series matrix has no expression table")
    expr = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t", index_col=0, quotechar='"')
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce")

    samples = pd.DataFrame(meta)
    if "geo_accession" in samples:
        samples = samples.set_index("geo_accession")
    samples = samples.reindex(expr.columns)
    return expr, samples, platform


def _symbol_from_gene_assignment(value: str) -> str:
    """HuGene/HuEx 'NM_001 // SYMBOL // desc /// ...' -> SYMBOL."""
    parts = [p.strip() for p in str(value).split("//")]
    return parts[1] if len(parts) > 1 and parts[1] not in ("", "---") else ""


def parse_platform_table(text: str) -> pd.Series:
    """Parse a GEO platform table (annot file or acc.cgi view=data) into probe -> symbol."""
    lines = [l for l in text.splitlines() if l and l[0] not in "^!#"]
    table = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", dtype=str, on_bad_lines="skip")
    table = table.set_index(table.columns[0])
    table.index = table.index.astype(str)
    lower = {c.lower(): c for c in table.columns}
    for name in ("gene symbol", "gene_symbol", "symbol", "identifier"):
        if name in lower:
            sym = table[lower[name]].fillna("")
            break
    else:
        if "gene_assignment" in lower:
            sym = table[lower["gene_assignment"]].fillna("").map(_symbol_from_gene_assignment)
        else:
            raise ValueError(f"no gene symbol column among {list(table.columns)[:15]}")
    sym = sym.astype(str).str.strip().replace({"---": ""})
    return sym


def platform_symbols(gpl: str, cache_dir: Path = DEFAULT_CACHE) -> pd.Series:
    """Probe -> gene symbol for a GEO platform. Tries the curated annot file first."""
    annot = f"{GEO_FTP}/platforms/{_prefix_dir(gpl)}/{gpl}/annot/{gpl}.annot.gz"
    try:
        return parse_platform_table(fetch(annot, cache_dir).decode("utf-8", "replace"))
    except requests.HTTPError:
        url = f"{GEO_ACC}?acc={gpl}&targ=self&form=text&view=data"
        return parse_platform_table(fetch(url, cache_dir).decode("utf-8", "replace"))


# --------------------------------------------------------------- processing

def needs_log2(expr: pd.DataFrame) -> bool:
    """GEO2R's heuristic for whether values are on a linear scale."""
    qx = np.nanquantile(expr.to_numpy(dtype=float), [0, 0.25, 0.5, 0.75, 0.99, 1.0])
    return bool(qx[4] > 100 or (qx[5] - qx[0] > 50 and qx[1] > 0))


def collapse_to_genes(expr: pd.DataFrame, symbols: pd.Series) -> pd.DataFrame:
    """One row per gene symbol, keeping the probe with the highest mean.

    Probes with no symbol or mapping to several genes ('A /// B') are dropped.
    """
    sym = symbols.reindex(expr.index).fillna("")
    keep = (sym != "") & ~sym.str.contains("///")
    e = expr[keep].copy()
    e["_symbol"] = sym[keep]
    e["_mean"] = e.drop(columns="_symbol").mean(axis=1)
    e = e.sort_values("_mean", ascending=False).drop_duplicates("_symbol")
    return e.set_index("_symbol").drop(columns="_mean").rename_axis("gene").sort_index()


def sample_text(samples: pd.DataFrame) -> pd.Series:
    """All descriptive metadata for each sample joined into one searchable string."""
    cols = [c for c in samples.columns
            if c.startswith(("title", "source_name", "characteristics", "description"))]
    return samples[cols].fillna("").astype(str).agg(" | ".join, axis=1)


def assign_groups(samples: pd.DataFrame, rules: dict[str, str]) -> pd.Series:
    """Label each sample with the first group whose regex matches its metadata.

    `rules` maps group name -> case-insensitive regex, e.g.
    {"case": r"alzheimer", "control": r"normal|control"}. Samples matching no
    rule (or matching more than one) get NaN and are left out of comparisons.
    """
    text = sample_text(samples)
    hits = pd.DataFrame({g: text.str.contains(rx, case=False, regex=True) for g, rx in rules.items()})
    group = pd.Series(np.nan, index=samples.index, dtype=object)
    single = hits.sum(axis=1) == 1
    group[single] = hits[single].idxmax(axis=1)
    return group


def load_table(expression: str | Path, samples: str | Path, name: str | None = None) -> ExpressionDataset:
    """Load your own data: an expression table and a sample sheet.

    expression: CSV/TSV, first column gene symbols, one column per sample
        (log2 or linear; linear values are log2-transformed). Repeated gene
        symbols are collapsed to the row with the highest mean.
    samples: CSV/TSV with a column of sample names matching the expression
        columns ("sample", or the first column) and a "group" column.
    """
    sep = lambda p: "," if str(p).lower().endswith((".csv", ".csv.gz")) else "\t"
    expr = pd.read_csv(expression, sep=sep(expression), index_col=0)
    expr.index = expr.index.astype(str).str.strip()
    expr = expr.apply(pd.to_numeric, errors="coerce")
    sheet = pd.read_csv(samples, sep=sep(samples), dtype=str)
    sheet = sheet.set_index("sample" if "sample" in sheet.columns else sheet.columns[0])
    if "group" not in sheet.columns:
        raise ValueError(f"{samples} needs a 'group' column; found {list(sheet.columns)}")
    missing = [s for s in sheet.index if s not in expr.columns]
    if missing:
        raise ValueError(f"samples not found in {expression}: {missing[:5]}")
    expr = expr[list(sheet.index)]
    notes = [f"source: {expression}", f"rows: {len(expr)}"]
    if needs_log2(expr):
        expr = np.log2(expr.clip(lower=1))
        notes.append("values looked linear; applied log2")
    symbols = pd.Series(expr.index.to_numpy(), index=[f"row{i}" for i in range(len(expr))])
    expr.index = symbols.index  # gene symbols may repeat; collapse by row position
    genes = collapse_to_genes(expr, symbols)
    notes.append(f"genes: {len(genes)}")
    notes.append(f"group counts: {sheet['group'].value_counts().to_dict()}")
    return ExpressionDataset(name or Path(expression).name.split(".")[0], "user", genes, sheet, notes)


def load_series(
    gse: str,
    rules: dict[str, str],
    platform: str | None = None,
    sample_filter: str | None = None,
    exclude: str | None = None,
    cache_dir: Path = DEFAULT_CACHE,
) -> ExpressionDataset:
    """Download one GEO series and return a log2 gene x sample matrix with group labels.

    platform: for multi-platform series, which GPL to use.
    sample_filter: optional regex; only samples whose metadata matches are kept
        (e.g. one brain region of a multi-region study).
    exclude: optional regex; samples whose metadata matches are dropped
        (e.g. cell lines mixed in with tissue samples).
    """
    urls = series_matrix_urls(gse, cache_dir)
    if platform:
        urls = [u for u in urls if f"-{platform}_" in u] or urls
    if len(urls) > 1:
        raise ValueError(f"{gse} has several platforms; pass platform= one of {urls}")
    expr, samples, gpl = parse_series_matrix(fetch(urls[0], cache_dir).decode("utf-8", "replace"))
    notes = [f"source: {urls[0]}", f"probes: {len(expr)}"]

    text = sample_text(samples)
    keep = pd.Series(True, index=samples.index)
    if sample_filter:
        keep &= text.str.contains(sample_filter, case=False, regex=True)
    if exclude:
        keep &= ~text.str.contains(exclude, case=False, regex=True)
    if not keep.all():
        samples, expr = samples[keep], expr.loc[:, keep.to_numpy()]
        notes.append(f"filters kept {int(keep.sum())} of {len(keep)} samples")

    if needs_log2(expr):
        expr = np.log2(expr.clip(lower=1))
        notes.append("values looked linear; applied log2")

    genes = collapse_to_genes(expr, platform_symbols(gpl, cache_dir))
    notes.append(f"genes after collapsing probes: {len(genes)}")

    samples = samples.copy()
    samples["group"] = assign_groups(samples, rules)
    counts = samples["group"].value_counts(dropna=False).to_dict()
    notes.append(f"group counts: {counts}")
    return ExpressionDataset(gse, gpl or platform or "", genes, samples, notes)
