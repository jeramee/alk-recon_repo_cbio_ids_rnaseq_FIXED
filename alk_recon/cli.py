from __future__ import annotations

"""CLI entrypoint.

Unit tests call:
  python -m alk_recon.cli --input <variants.tsv> --out <outdir>

So we support --out (preferred) and --outdir (legacy).
"""

import argparse
from pathlib import Path

from alk_recon.pipeline import run_pipeline


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="alk_recon")
    ap.add_argument("--input", required=True, help="Input variants TSV")
    ap.add_argument("--out", dest="out", default=None, help="Output directory")
    ap.add_argument("--outdir", dest="outdir", default=None, help="Legacy output directory")

    args = ap.parse_args(argv)

    out = args.out or args.outdir or "out"
    run_pipeline(input_path=Path(args.input), outdir=Path(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
