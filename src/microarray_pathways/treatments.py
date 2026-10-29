"""Treatment hypotheses from the pathway results.

Two independent routes, both producing ranked lists of existing drugs to look at,
not recommendations:

  targets    drug-target mapping. Take the enzymes driving each significantly
             shifted pathway (GSEA leading edge, gene itself changed in the same
             direction), look up drugs that act on them in KEGG DRUG, and keep the
             ones whose action would push the gene back (inhibitors of enzymes that
             went up, activators of enzymes that went down).

  reversal   signature reversal (the Connectivity Map idea). The LINCS L1000
             consensus signatures list, for ~5,400 compounds, the genes each one
             pushes up and down in cultured human cells. A compound whose "up"
             genes are low in the case group and whose "down" genes are high
             would, if it acted the same way in the patient tissue, push the
             expression profile back toward control. Scored on all genes and on
             the genes of the shifted metabolic pathways only.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from .geo import DEFAULT_CACHE, fetch
from .pathways import KEGG_REST, parse_kegg_symbols

ENRICHR_LIBRARY = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}"
L1000_LIBRARY = "LINCS_L1000_Chem_Pert_Consensus_Sigs"

INHIBIT = re.compile(r"inhibitor|antagonist|blocker|inverse agonist|antimetabolite|inactivator|suppress", re.I)
# "stimulant" and "inducer" describe effects, not what the drug does to its target, so they are left out
ACTIVATE = re.compile(r"(?<!inverse )agonist|activator|potentiator|opener", re.I)


# ------------------------------------------------------------ shared helpers

def pathway_genes(gsea: pd.DataFrame, de: pd.DataFrame, fdr: float = 0.05,
                  gene_fdr: float = 0.05) -> pd.DataFrame:
    """Genes that drive the significantly shifted pathways.

    A gene qualifies when it is in the leading edge of a pathway with GSEA FDR < `fdr`
    and is itself changed (FDR < `gene_fdr`) in the same direction as the pathway.
    Returns one row per gene: log2fc, t, fdr, direction and the pathways it drives.
    """
    rows = defaultdict(list)
    for r in gsea[gsea["fdr"] < fdr].itertuples():
        sign = 1 if r.direction == "up" else -1
        for g in str(r.leading_edge).split(";"):
            if g in de.index and de.at[g, "fdr"] < gene_fdr and np.sign(de.at[g, "t"]) == sign:
                rows[g].append(r.pathway)
    out = de.loc[list(rows), ["log2fc", "t", "fdr"]].copy()
    out["direction"] = np.where(out["t"] > 0, "up", "down")
    out["pathways"] = [rows[g] for g in out.index]
    return out.sort_values("t", key=abs, ascending=False)


def _drug_key(name: str) -> str:
    """Normalize a drug name for matching across databases: 'Metformin hydrochloride (JP19)' -> 'metformin'."""
    name = re.sub(r"\(.*?\)", "", name).strip().lower()
    name = re.sub(r"\b(" + _SALTS + r")\b", "", name)
    return re.sub(r"[^a-z0-9]", "", name)


_SALTS = (r"(?:di|tri|mono)?(?:hydrochloride|hydrobromide|hydrate|hemihydrate|hemipentahydrate|heptahydrate|"
          r"sesquihydrate|sodium|disodium|potassium|dipotassium|calcium|magnesium|mesylate|maleate|sulfate|"
          r"acetate|citrate|tartrate|phosphate|besylate|tosylate|fumarate|succinate|bromide|chloride|"
          r"anhydrous|hcl|tromethamine|salt|meglumine|lysine)")


# ------------------------------------------------------------ drug targets (KEGG DRUG)

def parse_kegg_drug(entry: str) -> dict:
    """Parse one KEGG DRUG flat-file entry into the fields used here."""
    fields: dict[str, list[str]] = defaultdict(list)
    key = None
    for line in entry.splitlines():
        if line.startswith("///"):
            break
        if line[:12].strip():
            key = line[:12].strip()
        if key:
            fields[key].append(line[12:].strip())
    names = [n.rstrip(";").strip() for n in fields.get("NAME", []) if n.strip()]
    entry_id = (fields.get("ENTRY") or [""])[0].split()[0] if fields.get("ENTRY") else ""
    targets = []
    for t in fields.get("TARGET", []):
        m = re.match(r"([A-Za-z0-9\-_.]+)\s*\[HSA:([\d ]+)\]", t)
        if m:
            targets.append({"symbol": m.group(1), "hsa": [f"hsa:{x}" for x in m.group(2).split()]})
    remark = " ".join(fields.get("REMARK", []))
    efficacy = " ".join(fields.get("EFFICACY", []))
    return {
        "id": entry_id,
        "name": re.sub(r"\s*\((?:[A-Z/0-9]+|TN)\)\s*$", "", names[0]) if names else entry_id,
        "names": names,
        "efficacy": efficacy,
        "classes": " ".join(fields.get("CLASS", [])),
        "targets": targets,
        "marketed": any("(TN)" in n for n in names) or "ATC code" in remark,
        "atc": " ".join(re.findall(r"ATC code: ([A-Z0-9 ]+)", remark)).strip(),
    }


def drug_action(drug: dict) -> str:
    """'inhibits', 'activates' or 'unclear', from KEGG's efficacy and class text."""
    text = f"{drug['efficacy']} {drug['classes']}"
    inh, act = bool(INHIBIT.search(text)), bool(ACTIVATE.search(text))
    if inh and not act:
        return "inhibits"
    if act and not inh:
        return "activates"
    return "unclear"


def kegg_gene_drugs(cache_dir: Path = DEFAULT_CACHE) -> dict[str, set[str]]:
    """Gene symbol -> KEGG DRUG ids that act on it (KEGG /link/drug/hsa)."""
    symbol = parse_kegg_symbols(fetch(f"{KEGG_REST}/list/hsa", cache_dir).decode())
    links = defaultdict(set)
    for line in fetch(f"{KEGG_REST}/link/drug/hsa", cache_dir).decode().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] in symbol:
            links[symbol[parts[0]]].add(parts[1].replace("dr:", ""))
    return dict(links)


def kegg_drug_entries(ids: set[str], cache_dir: Path = DEFAULT_CACHE, batch: int = 10) -> dict[str, dict]:
    """Parsed KEGG DRUG entries for `ids`, cached in one JSON file (KEGG serves 10 per request)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "kegg_drug_entries.json"
    cache = json.loads(path.read_text()) if path.exists() else {}
    todo = sorted(set(ids) - set(cache))
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        resp = requests.get(f"{KEGG_REST}/get/" + "+".join(f"dr:{d}" for d in chunk), timeout=120)
        resp.raise_for_status()
        for text in resp.text.split("///"):
            if text.strip():
                d = parse_kegg_drug(text.strip() + "\n///")
                cache[d["id"]] = d
        for d in chunk:  # remember ids KEGG did not return so they are not re-requested
            cache.setdefault(d, None)
        path.write_text(json.dumps(cache))
    return {d: cache[d] for d in ids if cache.get(d)}


def target_candidates(genes: pd.DataFrame, gene_drugs: dict[str, set[str]],
                      drugs: dict[str, dict]) -> pd.DataFrame:
    """Rank drugs acting on pathway-driving genes.

    `genes` comes from pathway_genes(). A drug "opposes" the change if it inhibits a
    gene that went up or activates one that went down. Ranked by: opposes, then
    marketed, then the summed |t| of the genes it hits.
    """
    hits = defaultdict(list)
    for g in genes.index:
        for d in gene_drugs.get(g, ()):
            if d in drugs:
                hits[d].append(g)
    rows = []
    for d, gs in hits.items():
        drug = drugs[d]
        action = drug_action(drug)
        dirs = {genes.at[g, "direction"] for g in gs}
        opposes = (action == "inhibits" and "up" in dirs) or (action == "activates" and "down" in dirs)
        paths = sorted({p for g in gs for p in genes.at[g, "pathways"]})
        rows.append({
            "drug": drug["name"], "kegg_id": d, "action": action,
            "opposes_change": opposes, "marketed": drug["marketed"],
            "targets": ", ".join(f"{g} ({genes.at[g, 'log2fc']:+.2f})" for g in gs),
            "n_targets": len(gs),
            "score": float(sum(abs(genes.at[g, "t"]) for g in gs)),
            "pathways": "; ".join(paths), "efficacy": drug["efficacy"], "atc": drug["atc"],
        })
    cols = ["drug", "kegg_id", "action", "opposes_change", "marketed", "targets", "n_targets", "score",
            "pathways", "efficacy", "atc"]
    res = pd.DataFrame(rows, columns=cols)
    if res.empty:
        return res
    res = res.sort_values(["opposes_change", "marketed", "score"], ascending=False)
    # one row per active ingredient: KEGG lists each salt and hydrate separately
    res["drug"] = res["drug"].map(lambda n: re.sub(r"\s+" + _SALTS + r"\b.*$", "", n, flags=re.I))
    key = res["drug"].map(_drug_key)
    ids = res.groupby(key, sort=False)["kegg_id"].agg(" ".join)
    res = res[~key.duplicated()].copy()
    res["kegg_id"] = res["drug"].map(_drug_key).map(ids)
    return res.reset_index(drop=True)


# ------------------------------------------------------------ signature reversal (LINCS L1000)

def parse_signature_library(gmt: str) -> dict[str, dict[str, set[str]]]:
    """Enrichr up/down GMT ('Afatinib Up', 'Afatinib Down') -> {compound: {'up': set, 'down': set}}."""
    sigs: dict[str, dict[str, set[str]]] = defaultdict(lambda: {"up": set(), "down": set()})
    for line in gmt.splitlines():
        row = line.split("\t")
        if len(row) < 3:
            continue
        m = re.fullmatch(r"(.+?)\s+(up|down)", row[0].strip(), re.I)
        if not m:
            continue
        sigs[m.group(1)][m.group(2).lower()] |= {g.split(",")[0] for g in row[2:] if g}
    return {k: v for k, v in sigs.items() if v["up"] and v["down"]}


def l1000_signatures(cache_dir: Path = DEFAULT_CACHE) -> dict[str, dict[str, set[str]]]:
    return parse_signature_library(fetch(ENRICHR_LIBRARY.format(L1000_LIBRARY), cache_dir).decode())


def _signature_matrix(sigs: dict[str, dict[str, set[str]]], genes: pd.Index, min_genes: int):
    """Sparse compounds x genes matrix with +1/n_up on up genes and -1/n_down on down genes."""
    from scipy import sparse

    pos = {g: i for i, g in enumerate(genes)}
    names, rows, cols, vals, n_up, n_down = [], [], [], [], [], []
    for name, sig in sigs.items():
        up = [pos[g] for g in sig["up"] if g in pos]
        dn = [pos[g] for g in sig["down"] if g in pos and g not in sig["up"]]
        if len(up) < min_genes or len(dn) < min_genes:
            continue
        r = len(names)
        names.append(name)
        n_up.append(len(up))
        n_down.append(len(dn))
        rows += [r] * (len(up) + len(dn))
        cols += up + dn
        vals += [1 / len(up)] * len(up) + [-1 / len(dn)] * len(dn)
    m = sparse.csr_matrix((vals, (rows, cols)), shape=(len(names), len(genes)))
    return names, m, n_up, n_down


def connectivity(t: pd.Series, sigs: dict[str, dict[str, set[str]]], universe: set[str] | None = None,
                 min_genes: int = 5) -> pd.DataFrame:
    """How well each compound signature matches the case vs control t statistics.

    score = mean t of the compound's up genes minus mean t of its down genes, over
    genes in `universe` (default: all measured), scaled to a z-score as if the
    genes were drawn at random. Negative z: the compound moves genes opposite to
    the disease (a reverser); positive: it mimics it.

    Genes are not independent (co-regulated programs such as the cell cycle move
    together), so z is a ranking score, not a calibrated test: compare compounds by
    `percentile` (share of all compounds that reverse less strongly) rather than by
    a p-value.
    """
    t = t.dropna()
    genes = t.index if universe is None else t.index[t.index.isin(universe)]
    names, m, n_up, n_down = _signature_matrix(sigs, genes, min_genes)
    cols = ["compound", "z", "percentile", "n_up", "n_down"]
    if not names:
        return pd.DataFrame(columns=cols)
    obs = m @ t[genes].to_numpy()
    z = obs / (t[genes].std() * np.sqrt(1 / np.array(n_up) + 1 / np.array(n_down)))
    res = pd.DataFrame({"compound": names, "z": z, "n_up": n_up, "n_down": n_down})
    res["percentile"] = 100 * (1 - (res["z"].rank(method="max") - 1) / len(res))
    return res[cols]


def reversal_candidates(de: pd.DataFrame, sigs: dict[str, dict[str, set[str]]], metabolic: set[str],
                        driver_genes: pd.DataFrame | None = None, clinical: set[str] = frozenset(),
                        min_genes: int = 5) -> pd.DataFrame:
    """Compounds whose L1000 signature opposes the disease, on metabolic genes and on all genes.

    `metabolic` is every gene in any metabolic gene set. Ranked by the metabolic z;
    the all-gene z shows whether the compound also reverses the wider profile.
    `reversed_genes` lists pathway-driving genes the compound pushes back (arrow =
    the way the compound moves it). `clinical` holds _drug_key() names of clinical
    drugs (kegg_clinical_names), used to flag compounds that are not just research tools.
    """
    met = connectivity(de["t"], sigs, universe=metabolic, min_genes=min_genes)
    allg = connectivity(de["t"], sigs, min_genes=min_genes)
    res = met.merge(allg[["compound", "z", "percentile"]], on="compound", how="left",
                    suffixes=("_metabolic", "_all"))
    res = res.rename(columns={"n_up": "n_up_metabolic", "n_down": "n_down_metabolic"})
    t = de["t"]
    focus = set(driver_genes.index) if driver_genes is not None and len(driver_genes) else metabolic

    def reversed_genes(name: str) -> str:
        sig = sigs[name]
        raised = [g for g in sig["up"] & focus if g in t.index and t[g] < 0]
        lowered = [g for g in sig["down"] & focus if g in t.index and t[g] > 0]
        picked = sorted(raised + lowered, key=lambda g: -abs(t[g]))
        return ", ".join(f"{g}{'↑' if g in raised else '↓'}" for g in picked[:12])

    res["clinical_drug"] = res["compound"].map(_drug_key).isin(clinical)
    res["reversed_genes"] = res["compound"].map(reversed_genes) if len(res) else []
    cols = ["compound", "clinical_drug", "z_metabolic", "percentile_metabolic", "z_all", "percentile_all",
            "n_up_metabolic", "n_down_metabolic", "reversed_genes"]
    return res[cols].sort_values("z_metabolic").reset_index(drop=True)


def kegg_clinical_names(cache_dir: Path = DEFAULT_CACHE) -> set[str]:
    """_drug_key() of every name of KEGG DRUG entries with an official drug name (INN, USAN, JAN),
    a pharmacopoeia listing or a trade name, i.e. compounds used or tested in people."""
    keys = set()
    for line in fetch(f"{KEGG_REST}/list/drug", cache_dir).decode().splitlines():
        parts = line.split("\t", 1)
        if len(parts) < 2 or not re.search(r"\((?:TN|[A-Z/]*(?:JP|USP|INN|USAN)[A-Z/0-9]*)\)", parts[1]):
            continue
        keys.update(k for k in (_drug_key(n) for n in parts[1].split(";")) if len(k) > 3)
    return keys


# ------------------------------------------------------------ both

def treatment_hypotheses(de: pd.DataFrame, gsea: pd.DataFrame, gene_sets: dict[str, set[str]],
                         out_dir: str | Path,
                         cache_dir: Path = DEFAULT_CACHE, fdr: float = 0.05) -> dict[str, pd.DataFrame]:
    """Run both routes, write drugs_targets.tsv, drugs_reversal.tsv and drug_target_genes.tsv."""
    out = Path(out_dir)
    genes = pathway_genes(gsea, de, fdr=fdr)

    gene_drugs = kegg_gene_drugs(cache_dir)
    wanted = {d for g in genes.index for d in gene_drugs.get(g, ())}
    targets = target_candidates(genes, gene_drugs, kegg_drug_entries(wanted, cache_dir))

    metabolic = set().union(*gene_sets.values()) if gene_sets else set()
    reversal = reversal_candidates(de, l1000_signatures(cache_dir), metabolic, genes,
                                   kegg_clinical_names(cache_dir))

    # cross-reference: compounds both routes pick up
    rev_z = {_drug_key(c): z for c, z in zip(reversal["compound"], reversal["z_metabolic"])}
    targets["l1000_z_metabolic"] = [rev_z.get(_drug_key(n), np.nan) for n in targets["drug"]]
    opposing = set(targets.loc[targets["opposes_change"].astype(bool), "drug"].map(_drug_key))
    reversal["known_target_drug"] = reversal["compound"].map(_drug_key).isin(opposing)

    genes.assign(pathways=genes["pathways"].map("; ".join)).to_csv(
        out / "drug_target_genes.tsv", sep="\t", float_format="%.4g", index_label="gene")
    targets.to_csv(out / "drugs_targets.tsv", sep="\t", index=False, float_format="%.4g")
    reversal.to_csv(out / "drugs_reversal.tsv", sep="\t", index=False, float_format="%.4g")
    return {"genes": genes, "targets": targets, "reversal": reversal}
