"""Hadronic-vacuum-polarisation plot.

Left: log-log plot of the renormalised HVP Pi(p^2) - Pi(0) (one colour, unit
charge). Right: [Pi(p^2) - Pi(0)]/p^2, whose p^2 -> 0 limit is the slope
Pi'(0); this is the panel to judge low-energy resolution.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .reader import Run  # noqa: E402


def plot_hvp(run: Run, outdir: str) -> str | None:
    if not run.has_hvp:
        return None
    os.makedirs(outdir, exist_ok=True)

    p_sq = run.hvp_p_sq
    pi = run.hvp_pi
    mask = pi > 0.0

    fig, (ax_log, ax_slope) = plt.subplots(1, 2, figsize=(11, 4))
    fig.suptitle(f"HVP  —  {run.label}")

    ax_log.loglog(p_sq[mask], pi[mask], "o-", color="black", ms=4)
    ax_log.set_xlabel(r"$p^2\ [\mathrm{GeV}^2]$")
    ax_log.set_ylabel(r"$\hat\Pi(p^2) = \Pi(p^2) - \Pi(0)$")

    ax_slope.semilogx(p_sq, pi / p_sq, "o-", color="black", ms=4)
    ax_slope.set_xlabel(r"$p^2\ [\mathrm{GeV}^2]$")
    ax_slope.set_ylabel(r"$\hat\Pi(p^2)\,/\,p^2\ [\mathrm{GeV}^{-2}]$")

    path = os.path.join(outdir, "hvp.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path
