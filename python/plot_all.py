#!/usr/bin/env python3
"""Generate every plot for a Quark_Photon_Vertex HDF5 run.

Usage:
    python -m plot_all output_charm_dse_pv.h5 [--outdir plots/]
    python plot_all.py   output_charm_dse_pv.h5 --outdir plots/charm/

Produces the f/g surfaces (fg_*.pdf), the four WTI comparisons (wti1..4.pdf),
and the HVP log-log plot (hvp.pdf, if the run contains /hvp).
"""

from __future__ import annotations

import argparse
import os
import sys

# Allow running both as a module (python -m plot_all) and as a script.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qpv import open_run
from qpv.plot_fg import plot_all_fg
from qpv.plot_wti import plot_all_wti
from qpv.plot_hvp import plot_hvp


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot a Quark_Photon_Vertex HDF5 run.")
    ap.add_argument("h5file", help="path to output_<flavour>[_dse][_pv].h5")
    ap.add_argument("--outdir", default=None,
                    help="output directory (default: plots_<run tag>/)")
    args = ap.parse_args()

    run = open_run(args.h5file)
    outdir = args.outdir or ("plots_" + os.path.splitext(os.path.basename(args.h5file))[0])
    os.makedirs(outdir, exist_ok=True)

    print(f"Run: {run.label}  ->  {outdir}/")
    fg = plot_all_fg(run, outdir)
    print(f"  {len(fg)} f/g surfaces")
    wti = plot_all_wti(run, outdir)
    print(f"  {len(wti)} WTI plots")
    hvp = plot_hvp(run, outdir)
    print("  hvp.pdf" if hvp else "  (no /hvp group — skipped HVP plot)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
