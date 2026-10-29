"""Metabolic pathway gene sets from KEGG and Reactome.

KEGG: human pathways from the KEGG REST API, kept if they belong to the
"Metabolism" category (map numbers below 01100). The global overview maps
(01100 "Metabolic pathways", 01200 "Carbon metabolism", ...) are left out
because they span hundreds of genes and swamp the specific pathways.

Reactome: the human pathway GMT, kept if the pathway sits under the top-level
"Metabolism" event (R-HSA-1430728) in the Reactome hierarchy.

Gene sets are keyed by a readable name, e.g. "KEGG: Glycolysis / Gluconeogenesis",
and map to sets of HGNC gene symbols (matching the loader's gene index).
"""

from __future__ import annotations

import io
import re
import zipfile
from collections import defaultdict
from pathlib import Path

from .geo import DEFAULT_CACHE, fetch

KEGG_REST = "https://rest.kegg.jp"
REACTOME_DL = "https://reactome.org/download/current"
REACTOME_METABOLISM = "R-HSA-1430728"


def _tsv(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line.strip()]


def parse_kegg(pathway_list: str, links: str, genes: str, metabolic_only: bool = True) -> dict[str, set[str]]:
    """Build KEGG gene sets from the text of /list/pathway/hsa, /link/hsa/pathway and /list/hsa."""
    names = {}
    for pid, name, *_ in _tsv(pathway_list):
        pid = pid.replace("path:", "")
        names[pid] = re.sub(r" - Homo sapiens \(human\)$", "", name)

    symbol = parse_kegg_symbols(genes)
    sets: dict[str, set[str]] = defaultdict(set)
    for pid, gid, *_ in _tsv(links):
        pid = pid.replace("path:", "")
        if metabolic_only and not is_kegg_metabolic(pid):
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


def parse_reactome(gmt: str, relations: str | None = None, root: str = REACTOME_METABOLISM) -> dict[str, set[str]]:
    """Build Reactome gene sets from the GMT; if `relations` is given keep only pathways under `root`."""
    keep = None
    if relations is not None:
        children = defaultdict(list)
        for parent, child, *_ in _tsv(relations):
            children[parent].append(child)
        keep, stack = set(), [root]
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


def kegg_metabolic(cache_dir: Path = DEFAULT_CACHE) -> dict[str, set[str]]:
    get = lambda path: fetch(f"{KEGG_REST}/{path}", cache_dir).decode()
    return parse_kegg(get("list/pathway/hsa"), get("link/hsa/pathway"), get("list/hsa"))


def reactome_metabolic(cache_dir: Path = DEFAULT_CACHE) -> dict[str, set[str]]:
    raw = fetch(f"{REACTOME_DL}/ReactomePathways.gmt.zip", cache_dir)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        gmt = z.read(next(n for n in z.namelist() if n.endswith(".gmt"))).decode()
    rel = fetch(f"{REACTOME_DL}/ReactomePathwaysRelation.txt", cache_dir).decode()
    return parse_reactome(gmt, rel)


# Ribosomal proteins, proteasome subunits, nucleoporins and translation factors.
# Several Reactome "Metabolism" pathways (e.g. Selenocysteine synthesis, Regulation
# of ODC, IP transport) are mostly these, so they track ribosome/proteasome/nuclear
# pore expression rather than metabolic enzymes.
MACHINERY = re.compile(r"^(RP[LS]\d|RPLP|PSM[ABCDEFG]|NUP\d|SEC13|EEF|EIF)")


def machinery_fraction(genes: set[str]) -> float:
    return sum(bool(MACHINERY.match(g)) for g in genes) / len(genes) if genes else 0.0


def metabolic_gene_sets(sources: tuple[str, ...] = ("kegg", "reactome"),
                        cache_dir: Path = DEFAULT_CACHE,
                        max_machinery: float | None = 0.5) -> dict[str, set[str]]:
    """Metabolic gene sets; sets over `max_machinery` housekeeping genes are dropped (None keeps all)."""
    sets: dict[str, set[str]] = {}
    if "kegg" in sources:
        sets.update(kegg_metabolic(cache_dir))
    if "reactome" in sources:
        sets.update(reactome_metabolic(cache_dir))
    if max_machinery is not None:
        sets = {k: v for k, v in sets.items() if machinery_fraction(v) <= max_machinery}
    return sets
