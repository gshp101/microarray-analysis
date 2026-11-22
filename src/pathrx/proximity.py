"""Network proximity of drugs to the disease module (Guney et al. 2016, Nat Commun 7:10331).

The disease module is the set of genes driving the shifted metabolic pathways
(treatments.pathway_genes). For each drug, the closest distance d_c is the mean,
over its targets, of the shortest-path distance in the human protein interaction
network to the nearest module gene. A small d_c alone means little, since drugs
that hit well connected proteins sit close to everything, so d_c is compared with
reference distances between random protein sets of the same size and degree
(both the target set and the module are re-drawn) and reported as a z-score.
Negative z: the drug's targets sit closer to the module than chance.

Network: the human interactome curated by Menche et al. (2015, Science), the one
used in the Guney paper (~13,400 proteins, ~141,000 interactions, Entrez ids).
Drug targets: KEGG DRUG, as in the drug-target route.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csgraph

from .geo import DEFAULT_CACHE, fetch
from .pathways import KEGG_REST, parse_kegg_symbols
from .treatments import _SALTS, _drug_key, pathway_genes

INTERACTOME_URL = "https://raw.githubusercontent.com/emreg00/proximity/master/data/network/network.sif"
CLINICAL_NAME = re.compile(r"\((?:TN|[A-Z/]*(?:JP|USP|INN|USAN)[A-Z/0-9]*)\)")


class Network:
    """Undirected graph on gene symbols, reduced to its largest connected component."""

    def __init__(self, edges: list[tuple[str, str]]):
        nodes = sorted({g for e in edges for g in e})
        idx = {g: i for i, g in enumerate(nodes)}
        a = np.array([idx[x] for x, y in edges if x != y])
        b = np.array([idx[y] for x, y in edges if x != y])
        adj = sparse.coo_matrix((np.ones(len(a)), (a, b)), shape=(len(nodes), len(nodes))).tocsr()
        adj = ((adj + adj.T) > 0).astype(np.float32)
        _, labels = csgraph.connected_components(adj, directed=False)
        keep = np.flatnonzero(labels == np.bincount(labels).argmax())
        self.adj = adj[keep][:, keep].tocsr()
        self.nodes = [nodes[i] for i in keep]
        self.index = {g: i for i, g in enumerate(self.nodes)}
        self.degree = np.asarray(self.adj.sum(axis=1)).ravel().astype(int)

    def __len__(self) -> int:
        return len(self.nodes)

    def ids(self, genes) -> np.ndarray:
        return np.array(sorted({self.index[g] for g in genes if g in self.index}), dtype=int)

    def distance_to(self, sources: np.ndarray) -> np.ndarray:
        """Distance from every node to the nearest node of each source set.

        `sources` is a boolean (n_nodes, n_sets) matrix; one breadth-first search per
        column, run together as sparse matrix products. Returns int16 (n_nodes, n_sets).
        """
        sources = np.asarray(sources, dtype=bool)
        if sources.ndim == 1:
            sources = sources[:, None]
        dist = np.full(sources.shape, -1, dtype=np.int16)
        dist[sources] = 0
        frontier, step = sources.astype(np.float32), 0
        while frontier.any():
            step += 1
            reached = (self.adj @ frontier > 0) & (dist < 0)
            dist[reached] = step
            frontier = reached.astype(np.float32)
        return dist


def parse_sif(text: str, entrez_symbol: dict[str, str]) -> list[tuple[str, str]]:
    """'geneA <type> geneB' lines with Entrez ids -> symbol pairs (ids without a symbol are dropped)."""
    edges = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0] in entrez_symbol and parts[2] in entrez_symbol:
            edges.append((entrez_symbol[parts[0]], entrez_symbol[parts[2]]))
    return edges


def load_interactome(cache_dir: Path = DEFAULT_CACHE) -> Network:
    kegg = parse_kegg_symbols(fetch(f"{KEGG_REST}/list/hsa", cache_dir).decode())
    entrez = {k.split(":", 1)[1]: v for k, v in kegg.items()}  # KEGG human gene ids are Entrez ids
    return Network(parse_sif(fetch(INTERACTOME_URL, cache_dir).decode(), entrez))


def kegg_drug_targets(cache_dir: Path = DEFAULT_CACHE) -> pd.DataFrame:
    """Every KEGG DRUG entry with human targets, merged per active ingredient.

    Returns one row per drug: name, kegg_id (space separated), clinical (official drug
    name, pharmacopoeia or trade name) and targets (set of symbols).
    """
    symbol = parse_kegg_symbols(fetch(f"{KEGG_REST}/list/hsa", cache_dir).decode())
    targets = defaultdict(set)
    for line in fetch(f"{KEGG_REST}/link/drug/hsa", cache_dir).decode().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] in symbol:
            targets[parts[1].replace("dr:", "")].add(symbol[parts[0]])
    rows = []
    for line in fetch(f"{KEGG_REST}/list/drug", cache_dir).decode().splitlines():
        parts = line.split("\t", 1)
        did = parts[0].replace("dr:", "")
        if len(parts) < 2 or did not in targets:
            continue
        name = re.sub(r"\s*\((?:[A-Z/0-9]+|TN)\)\s*$", "", parts[1].split(";")[0].strip())
        name = re.sub(r"\s+" + _SALTS + r"\b.*$", "", name, flags=re.I)
        rows.append({"drug": name, "kegg_id": did,
                     "clinical": bool(CLINICAL_NAME.search(parts[1])), "targets": targets[did]})
    df = pd.DataFrame(rows, columns=["drug", "kegg_id", "clinical", "targets"])
    df["key"] = df["drug"].map(_drug_key)
    df = df[df["key"] != ""]
    return (df.groupby("key", sort=False)
              .agg(drug=("drug", lambda s: min(s, key=len)), kegg_id=("kegg_id", " ".join),
                   clinical=("clinical", "any"), targets=("targets", lambda s: set().union(*s)))
              .reset_index(drop=True))


def degree_bins(degree: np.ndarray, min_bin_size: int = 100) -> np.ndarray:
    """Bin label per node: nodes sorted by degree, each bin holds whole degree values and
    at least `min_bin_size` nodes (the last bin absorbs any remainder)."""
    values, counts = np.unique(degree, return_counts=True)
    label, size, by_value = 0, 0, {}
    for v, c in zip(values, counts):
        by_value[v] = label
        size += c
        if size >= min_bin_size:
            label, size = label + 1, 0
    if size and label:  # fold an undersized last bin into the previous one
        by_value = {v: min(b, label - 1) for v, b in by_value.items()}
    return np.array([by_value[d] for d in degree])


class DegreeSampler:
    """Draw random node sets matching a given set's degrees, bin by bin, without repeats."""

    def __init__(self, degree: np.ndarray, min_bin_size: int = 100):
        self.bin = degree_bins(degree, min_bin_size)
        self.members = [np.flatnonzero(self.bin == b) for b in range(self.bin.max() + 1)]

    def sample(self, nodes: np.ndarray, reps: int, rng: np.random.Generator) -> np.ndarray:
        """(reps, len(nodes)) node ids; column j is drawn from the degree bin of nodes[j]."""
        out = np.empty((reps, len(nodes)), dtype=int)
        bins = self.bin[nodes]
        for b in np.unique(bins):
            cols = np.flatnonzero(bins == b)
            pool = self.members[b]
            pick = np.argsort(rng.random((reps, len(pool))), axis=1)[:, :len(cols)]
            out[:, cols] = pool[pick]
        return out


def drug_proximity(net: Network, module: list[str], drugs: pd.DataFrame, reps: int = 1000,
                   min_bin_size: int = 100, seed: int = 0) -> pd.DataFrame:
    """Closest distance and z-score from each drug's targets to the module genes.

    `drugs` has drug, kegg_id, clinical and targets columns (kegg_drug_targets). Drugs
    with no target in the network are left out.
    """
    cols = ["drug", "kegg_id", "clinical", "z", "d_closest", "d_random_mean", "d_random_sd",
            "n_targets", "n_targets_in_module", "targets"]
    mod = net.ids(module)
    if len(mod) == 0 or drugs.empty:
        return pd.DataFrame(columns=cols)
    rng = np.random.default_rng(seed)
    sampler = DegreeSampler(net.degree, min_bin_size)

    mask = np.zeros(len(net), dtype=bool)
    mask[mod] = True
    d_mod = net.distance_to(mask)[:, 0].astype(float)

    # reference: distance to random degree-matched modules, one column per repetition
    rand_mod = sampler.sample(mod, reps, rng)
    masks = np.zeros((len(net), reps), dtype=bool)
    masks[rand_mod, np.arange(reps)[:, None]] = True
    d_rand = net.distance_to(masks).T  # (reps, nodes)

    module_set = set(module)
    rows, reference = [], {}
    for r in drugs.itertuples(index=False):
        t = net.ids(r.targets)
        if len(t) == 0:
            continue
        d = d_mod[t].mean()
        key = tuple(t)  # many drugs share a target set; draw its reference once
        if key not in reference:
            reference[key] = d_rand[np.arange(reps)[:, None], sampler.sample(t, reps, rng)].mean(axis=1)
        ref = reference[key]
        sd = ref.std()
        in_mod = sorted(set(r.targets) & module_set)
        rows.append({
            "drug": r.drug, "kegg_id": r.kegg_id, "clinical": r.clinical,
            "z": (d - ref.mean()) / sd if sd > 0 else np.nan,
            "d_closest": d, "d_random_mean": ref.mean(), "d_random_sd": sd,
            "n_targets": len(t), "n_targets_in_module": len(in_mod),
            "targets": ", ".join(in_mod + sorted(net.nodes[i] for i in t if net.nodes[i] not in module_set)),
        })
    return pd.DataFrame(rows, columns=cols).sort_values(["z", "d_closest"]).reset_index(drop=True)


def proximity_hypotheses(de: pd.DataFrame, gsea: pd.DataFrame, out_dir: str | Path,
                         genes: pd.DataFrame | None = None, targets: pd.DataFrame | None = None,
                         reversal: pd.DataFrame | None = None, cache_dir: Path = DEFAULT_CACHE,
                         fdr: float = 0.05, reps: int = 1000, seed: int = 0) -> dict:
    """Rank KEGG drugs by network proximity to the pathway-driving genes; write drugs_proximity.tsv.

    `genes`, `targets` and `reversal` are the treatment step's tables, if it ran: the
    module is taken from `genes`, and the other two are cross-referenced.
    """
    out = Path(out_dir)
    if genes is None:
        genes = pathway_genes(gsea, de, fdr=fdr)
    net = load_interactome(cache_dir)
    module = [g for g in genes.index if g in net.index]
    res = drug_proximity(net, module, kegg_drug_targets(cache_dir), reps=reps, seed=seed)

    if targets is not None and len(targets):
        opposing = set(targets.loc[targets["opposes_change"].astype(bool), "drug"].map(_drug_key))
        res["opposes_change"] = res["drug"].map(_drug_key).isin(opposing)
    if reversal is not None and len(reversal):
        rev_z = {_drug_key(c): z for c, z in zip(reversal["compound"], reversal["z_metabolic"])}
        res["l1000_z_metabolic"] = [rev_z.get(_drug_key(n), np.nan) for n in res["drug"]]
    res.to_csv(out / "drugs_proximity.tsv", sep="\t", index=False, float_format="%.4g")
    info = {"module_genes": len(module), "pathway_genes": len(genes), "network_nodes": len(net),
            "network_edges": int(net.adj.nnz // 2), "reps": reps}
    return {"proximity": res, "proximity_info": info}
