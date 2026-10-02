"""Hadronic-vacuum-polarisation plot.

Ports the old repo-root plot_hvp.py to the HDF5 layout: log-log plot of the
renormalised Pi(p^2) vs p^2, masking non-positive Pi (as the original did).
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

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.loglog(p_sq[mask], pi[mask], "o-", color="black", ms=4,
              label=r"$\tilde{\Pi}(p^2)$")
    ax.set_xlabel(r"$p^2\ [\mathrm{GeV}^2]$")
    ax.set_ylabel(r"$\tilde{\Pi}(p^2)$")
    ax.set_title(f"HVP  —  {run.label}")
    ax.legend(loc="best", frameon=False)

    path = os.path.join(outdir, "hvp.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path
