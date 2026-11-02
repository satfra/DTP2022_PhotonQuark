"""Quark-DSE propagator functions: dressing Z_q(p^2) and mass M_q(p^2).

Only present for runs solved with the quark DSE (-d). Both curves share one
axes: the IR plateau of M_q is the dynamically generated (constituent) mass, the
UV tail is the current mass.

Note: Z_q here is the code's vector dressing A(p^2). If you want the wave-function
renormalisation 1/A instead, invert `A` below.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .reader import Run  # noqa: E402


def plot_dse(run: Run, outdir: str) -> str | None:
    if not run.has_dse:
        return None
    os.makedirs(outdir, exist_ok=True)

    d = run.dse
    p_sq, Zq, Mq = d["p_sq"], d["A"], d["M"]
    order = np.argsort(p_sq)
    p_sq, Zq, Mq = p_sq[order], Zq[order], Mq[order]
    p = np.sqrt(p_sq)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(p, 1/Zq, color="black", label=r"$Z_q^{-1}$")
    ax.plot(p, Mq, color="red", label=r"$M_q$")
    ax.set_xlabel(r"$p\ [\mathrm{GeV}]$")
    ax.set_xscale("log")
    ax.set_ylim(0.0, 1.1)
    #ax.set_xlim(1e-2, 50)
    ax.set_title(f"Quark DSE — {run.label}")
    ax.legend(loc="best", frameon=False)

    path = os.path.join(outdir, "quark_dse.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path
