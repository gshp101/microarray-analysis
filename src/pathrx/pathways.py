"""Metabolic (and optionally signaling) pathway gene sets from KEGG and Reactome.

KEGG: human pathways from the KEGG REST API, kept if they belong to the
"Metabolism" category (map numbers below 01100). The global overview maps
(01100 "Metabolic pathways", 01200 "Carbon metabolism", ...) are left out
because they span hundreds of genes and swamp the specific pathways.

Reactome: the human pathway GMT, kept if the pathway sits under the top-level
"Metabolism" event (R-HSA-1430728) in the Reactome hierarchy.

Signaling sets, used with scope "signaling" or "all": KEGG pathways in the BRITE
"Signal transduction" group (MAPK, PI3K-Akt, VEGF, ...; br:br08901) and Reactome
pathways under "Signal Transduction" (R-HSA-162582). Kinase inhibitors and other
signaling drugs only reach the drug routes when these are included.

Gene sets are keyed by a readable name, e.g. "KEGG: Glycolysis / Gluconeogenesis",
and map to sets of HGNC gene symbols (matching the loader's gene index).
"""

from __future__ import annotations

import io
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Callable

from .geo import DEFAULT_CACHE, fetch

KEGG_REST = "https://rest.kegg.jp"
REACTOME_DL = "https://reactome.org/download/current"
REACTOME_METABOLISM = "R-HSA-1430728"
REACTOME_SIGNALING = "R-HSA-162582"
SCOPES = ("metabolic", "signaling", "all")


def _tsv(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line.strip()]


def parse_kegg(pathway_list: str, links: str, genes: str,
               keep: Callable[[str], bool] | None = None) -> dict[str, set[str]]:
    """Build KEGG gene sets from the text of /list/pathway/hsa, /link/hsa/pathway and /list/hsa.

    `keep(pathway_id)` picks the pathways (default: the Metabolism category)."""
    keep = keep or is_kegg_metabolic
    names = {}
    for pid, name, *_ in _tsv(pathway_list):
        pid = pid.replace("path:", "")
        names[pid] = re.sub(r" - Homo sapiens \(human\)$", "", name)

    symbol = parse_kegg_symbols(genes)
    sets: dict[str, set[str]] = defaultdict(set)
    for pid, gid, *_ in _tsv(links):
        pid = pid.replace("path:", "")
        if not keep(pid):
            continue
        if gid in symbol and pid in names:
            sets[f"KEGG: {names[pid]}"].add(symbol[gid])
    return dict(sets)


def parse_kegg_symbols(genes: str) -> dict[str, str]:
    """KEGG gene id (hsa:3101) -> official symbol, from the text of /list/hsa."""
    symbol = {}
    for row in _tsv(genes):
        # hsa:3101 <tab> CDS <tab> location <tab> "HK3, HKIII; hexokinase 3"  (older dumps omit type/location)
        # genes with no official symbol have only a description and no ";"
        gid, desc = row[0], row[-1]
        if ";" not in desc:
            continue
        sym = desc.split(";", 1)[0].split(",", 1)[0].strip()
        if sym:
            symbol[gid] = sym
    return symbol


def is_kegg_metabolic(pid: str) -> bool:
    m = re.fullmatch(r"[a-z]+(\d{5})", pid)
    return bool(m) and int(m.group(1)) < 1100


def parse_kegg_brite_group(brite: str, group: str = "Signal transduction") -> set[str]:
    """Five-digit map numbers listed under one B-level group of the KEGG pathway hierarchy (br08901)."""
    maps, inside = set(), False
    for line in brite.splitlines():
        if line.startswith("B"):
            inside = line[1:].strip() == group
        elif line.startswith("A"):
            inside = False
        elif inside and line.startswith("C"):
            m = re.match(r"C\s+(\d{5})", line)
            if m:
                maps.add(m.group(1))
    return maps


def parse_reactome(gmt: str, relations: str | None = None,
                   root: str | tuple[str, ...] = REACTOME_METABOLISM) -> dict[str, set[str]]:
    """Build Reactome gene sets from the GMT; if `relations` is given keep only pathways under `root`
    (one id or several)."""
    keep = None
    if relations is not None:
        children = defaultdict(list)
        for parent, child, *_ in _tsv(relations):
            children[parent].append(child)
        keep, stack = set(), [root] if isinstance(root, str) else list(root)
        while stack:
            node = stack.pop()
            if node not in keep:
                keep.add(node)
                stack.extend(children[node])

    sets = {}
    for row in _tsv(gmt):
        if len(row) < 3:
            continue
        name, rid, genes = row[0], row[1], row[2:]
        if not rid.startswith("R-HSA-") or (keep is not None and rid not in keep):
            continue
        sets[f"Reactome: {name}"] = {g for g in genes if g}
    return sets


def _kegg_number(pid: str) -> str:
    return re.sub(r"^[a-z]+", "", pid)


def kegg_sets(scope: str = "metabolic", cache_dir: Path = DEFAULT_CACHE) -> dict[str, set[str]]:
    get = lambda path: fetch(f"{KEGG_REST}/{path}", cache_dir).decode()
    signaling = parse_kegg_brite_group(get("get/br:br08901")) if scope != "metabolic" else set()
    keep = {"metabolic": is_kegg_metabolic,
            "signaling": lambda pid: _kegg_number(pid) in signaling,
            "all": lambda pid: is_kegg_metabolic(pid) or _kegg_number(pid) in signaling}[scope]
    return parse_kegg(get("list/pathway/hsa"), get("link/hsa/pathway"), get("list/hsa"), keep=keep)


def reactome_sets(scope: str = "metabolic", cache_dir: Path = DEFAULT_CACHE) -> dict[str, set[str]]:
    raw = fetch(f"{REACTOME_DL}/ReactomePathways.gmt.zip", cache_dir)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        gmt = z.read(next(n for n in z.namelist() if n.endswith(".gmt"))).decode()
    rel = fetch(f"{REACTOME_DL}/ReactomePathwaysRelation.txt", cache_dir).decode()
    roots = {"metabolic": (REACTOME_METABOLISM,), "signaling": (REACTOME_SIGNALING,),
             "all": (REACTOME_METABOLISM, REACTOME_SIGNALING)}[scope]
    return parse_reactome(gmt, rel, roots)


def kegg_metabolic(cache_dir: Path = DEFAULT_CACHE) -> dict[str, set[str]]:
    return kegg_sets("metabolic", cache_dir)


def reactome_metabolic(cache_dir: Path = DEFAULT_CACHE) -> dict[str, set[str]]:
    return reactome_sets("metabolic", cache_dir)


# Ribosomal proteins, proteasome subunits, nucleoporins and translation factors.
# Several Reactome "Metabolism" pathways (e.g. Selenocysteine synthesis, Regulation
# of ODC, IP transport) are mostly these, so they track ribosome/proteasome/nuclear
# pore expression rather than metabolic enzymes.
MACHINERY = re.compile(r"^(RP[LS]\d|RPLP|PSM[ABCDEFG]|NUP\d|SEC13|EEF|EIF)")


def machinery_fraction(genes: set[str]) -> float:
    return sum(bool(MACHINERY.match(g)) for g in genes) / len(genes) if genes else 0.0


def pathway_gene_sets(sources: tuple[str, ...] = ("kegg", "reactome"), scope: str = "metabolic",
                      cache_dir: Path = DEFAULT_CACHE,
                      max_machinery: float | None = 0.5) -> dict[str, set[str]]:
    """Gene sets for `scope` (metabolic, signaling or all); sets over `max_machinery` housekeeping
    genes are dropped (None keeps all)."""
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {SCOPES}, got {scope!r}")
    sets: dict[str, set[str]] = {}
    if "kegg" in sources:
        sets.update(kegg_sets(scope, cache_dir))
    if "reactome" in sources:
        sets.update(reactome_sets(scope, cache_dir))
    if max_machinery is not None:
        sets = {k: v for k, v in sets.items() if machinery_fraction(v) <= max_machinery}
    return sets


def metabolic_gene_sets(sources: tuple[str, ...] = ("kegg", "reactome"),
                        cache_dir: Path = DEFAULT_CACHE,
                        max_machinery: float | None = 0.5) -> dict[str, set[str]]:
    """Metabolic gene sets; sets over `max_machinery` housekeeping genes are dropped (None keeps all)."""
    return pathway_gene_sets(sources, "metabolic", cache_dir, max_machinery)
