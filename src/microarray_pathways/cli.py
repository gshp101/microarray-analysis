"""Command line entry point.

    python -m microarray_pathways list
    python -m microarray_pathways load alzheimers_brain --out data/
    python -m microarray_pathways load GSE12345 --group case='tumor' --group control='normal'
    python -m microarray_pathways report alzheimers_brain
    python -m microarray_pathways report GSE12345 --group tumor='tumou?r' --group normal='normal'
    python -m microarray_pathways report --expression expr.tsv --samples samples.tsv --case AD --control CTL
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import datasets, geo


def _cmd_list(_: argparse.Namespace) -> int:
    for spec in datasets.DATASETS.values():
        status = "verified" if spec.verified else "rules not yet verified"
        print(f"{spec.key:22s} {spec.accession:9s} {spec.case} vs {spec.control}  ({status})")
        print(f"{'':22s} {spec.description}")
    return 0


def _cmd_load(args: argparse.Namespace) -> int:
    try:
        spec = datasets.get(args.dataset)
    except KeyError:
        spec = None

    if args.group:
        rules = dict(g.split("=", 1) for g in args.group)
        accession, platform, sample_filter, exclude, expected = (
            args.dataset, args.platform, args.filter, args.exclude, {})
    elif spec:
        rules, accession = spec.rules, spec.accession
        platform = args.platform or spec.platform
        sample_filter = args.filter or spec.sample_filter
        exclude = args.exclude or spec.exclude
        expected = spec.expected
    else:
        print(f"{args.dataset} is not a curated dataset; pass --group NAME=REGEX twice", file=sys.stderr)
        return 2

    ds = geo.load_series(accession, rules, platform=platform, sample_filter=sample_filter,
                         exclude=exclude, cache_dir=Path(args.cache))
    out = ds.save(args.out)
    print("\n".join(ds.notes))

    counts = ds.group_counts().to_dict()
    mismatch = {g: (counts.get(g, 0), n) for g, n in expected.items() if counts.get(g, 0) != n}
    if mismatch:
        print("WARNING group sizes differ from the study (got, expected):", mismatch)
    print(f"saved to {out}")
    return 1 if mismatch else 0


def _resolve(args: argparse.Namespace):
    """Pick the dataset to analyze: own table, curated dataset, or any GSE with --group rules.

    Returns (dataset, case, control, label, source).
    """
    if args.expression:
        if not (args.samples and args.case and args.control):
            raise SystemExit("--expression needs --samples, --case and --control")
        ds = geo.load_table(args.expression, args.samples, name=args.name)
        return ds, args.case, args.control, args.name or ds.accession, f"user table {args.expression}"

    if not args.dataset:
        raise SystemExit("give a curated dataset, a GSE accession with --group, or --expression/--samples")
    try:
        spec = None if args.group else datasets.get(args.dataset)
    except KeyError:
        raise SystemExit(f"{args.dataset} is not a curated dataset; pass --group NAME=REGEX twice "
                         "(case group first)")
    if spec:
        accession, rules, case, control = spec.accession, spec.rules, spec.case, spec.control
        platform, sample_filter, exclude = spec.platform, spec.sample_filter, spec.exclude
    else:
        rules = dict(g.split("=", 1) for g in args.group)
        if len(rules) != 2:
            raise SystemExit("pass exactly two --group NAME=REGEX rules, case group first")
        accession = args.dataset.upper()
        case, control = list(rules)
        platform, sample_filter, exclude = None, None, None
    platform = args.platform or platform
    sample_filter = args.filter or sample_filter
    exclude = args.exclude or exclude

    saved = Path(args.data) / accession
    if spec and saved.exists() and not (args.platform or args.filter or args.exclude):
        ds = geo.ExpressionDataset.read(saved)
    else:
        ds = geo.load_series(accession, rules, platform=platform, sample_filter=sample_filter,
                             exclude=exclude, cache_dir=Path(args.cache))
        ds.save(args.data)
    source = f"GEO {accession}" + (f" ({spec.description})" if spec else "")
    return ds, case, control, f"{accession}: {case} vs {control}", source


def _cmd_report(args: argparse.Namespace) -> int:
    from . import analysis, pathways

    ds, case, control, label, source = _resolve(args)
    sets = pathways.metabolic_gene_sets(tuple(args.sources.split(",")), Path(args.cache))
    analysis.analyze(ds, case, control, sets, args.out, label=label,
                     permutations=args.permutations, source=source, drugs=not args.no_drugs,
                     cache_dir=Path(args.cache))
    out = Path(args.out) / ds.accession
    print((out / "summary.txt").read_text(), end="")
    print(f"report: {out / 'report.html'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="microarray_pathways")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="list curated datasets").set_defaults(func=_cmd_list)

    lp = sub.add_parser("load", help="download a GEO series and build a gene matrix")
    lp.add_argument("dataset", help="curated key (see list) or a GSE accession")
    lp.add_argument("--group", action="append", metavar="NAME=REGEX",
                    help="group rule, regex matched against sample metadata (repeat)")
    lp.add_argument("--platform", help="GPL id for multi-platform series")
    lp.add_argument("--filter", help="keep only samples whose metadata matches this regex")
    lp.add_argument("--exclude", help="drop samples whose metadata matches this regex")
    lp.add_argument("--out", default="data", help="output directory")
    lp.add_argument("--cache", default=str(geo.DEFAULT_CACHE), help="download cache directory")
    lp.set_defaults(func=_cmd_load)

    ap = sub.add_parser("report", aliases=["analyze"],
                        help="differential expression + metabolic pathways + HTML report")
    ap.add_argument("dataset", nargs="?", help="curated key (see list) or a GSE accession")
    ap.add_argument("--group", action="append", metavar="NAME=REGEX",
                    help="for any GSE: group rule matched against sample metadata; give two, case first")
    ap.add_argument("--platform", help="GPL id for multi-platform series")
    ap.add_argument("--filter", help="keep only samples whose metadata matches this regex")
    ap.add_argument("--exclude", help="drop samples whose metadata matches this regex")
    ap.add_argument("--expression", help="your own expression table (genes x samples, CSV/TSV)")
    ap.add_argument("--samples", help="sample sheet for --expression with 'sample' and 'group' columns")
    ap.add_argument("--case", help="case group name in the sample sheet")
    ap.add_argument("--control", help="control group name in the sample sheet")
    ap.add_argument("--name", help="name for the output folder and report title (own table)")
    ap.add_argument("--data", default="data", help="where loaded GEO datasets are saved")
    ap.add_argument("--out", default="results", help="output directory")
    ap.add_argument("--sources", default="kegg,reactome", help="pathway databases, comma separated")
    ap.add_argument("--permutations", type=int, default=1000, help="GSEA permutations")
    ap.add_argument("--no-drugs", action="store_true",
                    help="skip the treatment-hypothesis step (KEGG DRUG targets, LINCS L1000 reversal)")
    ap.add_argument("--cache", default=str(geo.DEFAULT_CACHE), help="download cache directory")
    ap.set_defaults(func=_cmd_report)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
