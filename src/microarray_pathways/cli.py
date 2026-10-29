"""Command line entry point.

    python -m microarray_pathways list
    python -m microarray_pathways load alzheimers_brain --out data/
    python -m microarray_pathways load GSE12345 --group case='tumor' --group control='normal'
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

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
