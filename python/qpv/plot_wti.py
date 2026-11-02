"""Ward-Takahashi-identity check — a single 2x2 grid figure (wti.pdf).

Reproduces the old plot_wtis.py layout on the HDF5 output. In the exact-WTI
limit the z-averaged g-functions match the propagator quantities:
    WTI1:  g1  ==  Sigma_A
    WTI2:  g2  ==  2 * Delta_A
    WTI3:  g3  == -2 * Delta_B
    WTI4:  g4  ==  0
Each panel overlays the two surfaces over (log10 Q^2, log10 k^2) and prints the
peak-normalised L-inf error max|g - target| / max|target| (or max|g4| for WTI4).

Signs/factors verified numerically against a converged run (least-squares fits
give g1 = +1.00 Sigma_A, g2 = +2 Delta_A, g3 = -2 Delta_B, all |corr| ~ 1). This
matches the original plot_wtis.py; the old gnuplot wti.gplt had "-g2" and a wrong
Delta_B data-file reference, both of which were bugs.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: E402,F401  (registers 3d projection)

from .reader import Run  # noqa: E402


def _plot_pair(ax, logQ2, logK2, Zg, Zt, g_label, t_label, wti_num) -> None:
    ax.plot_surface(logQ2, logK2, Zg, cmap="viridis", alpha=0.7,
                    linewidth=0, antialiased=False)
    if Zt is not None:
        ax.plot_surface(logQ2, logK2, Zt, cmap="plasma", alpha=0.7,
                        linewidth=0, antialiased=False)
    ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
    ax.set_ylabel(r"$k^2$ [GeV$^2$]")
    ax.set_zlabel(f"wti{wti_num}")
    ax.set_title(f"WTI {wti_num}: {g_label}" + (f" vs {t_label}" if t_label else ""))
    fmt = lambda v, _pos: f"$10^{{{int(round(v))}}}$"
    ax.xaxis.set_major_formatter(plt.FuncFormatter(fmt))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(fmt))

    handles = [Patch(facecolor=plt.cm.viridis(0.6), alpha=0.7, label=g_label)]
    if t_label:
        handles.append(Patch(facecolor=plt.cm.plasma(0.6), alpha=0.7, label=t_label))
    ax.legend(handles=handles, loc="upper left", fontsize=8)

    if Zt is not None:
        # Peak-normalised L_inf error: max|g - target| / max|target|. Pointwise
        # |g-t|/|t| blows up wherever the target crosses zero while g is small
        # but nonzero — invisible on the surface, but it inflates the number.
        peak = float(np.abs(Zt).max())
        max_rel = float(np.abs(Zg - Zt).max()) / peak if peak > 0 else float("nan")
        info = f"max rel. err. = {max_rel:.2e}"
    else:
        info = f"max |{g_label}| = {float(np.nanmax(np.abs(Zg))):.2e}"
    ax.text2D(0.95, 0.95, info, transform=ax.transAxes, fontsize=9,
              ha="right", va="top",
              bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))


def plot_wti(run: Run, outdir: str) -> str:
    """Write a single wti.pdf (2x2 grid of the four WTI checks)."""
    os.makedirs(outdir, exist_ok=True)

    # data surfaces, shape [q, k]
    g1, g2, g3, g4 = (np.real(run.fgz0(n)) for n in ("g1", "g2", "g3", "g4"))
    sigma_a = np.real(run.wz0("Sigma_A"))
    delta_a = np.real(run.wz0("Delta_A"))
    delta_b = np.real(run.wz0("Delta_B"))

    # matching [q, k] meshgrids of the log10 axes (indexing="ij")
    logQ2, logK2 = np.meshgrid(np.log10(run.q_grid), np.log10(run.k_sq), indexing="ij")

    fig = plt.figure(figsize=(14, 11))
    fig.suptitle(f"WTI checks — {run.label}", fontsize=14)

    _plot_pair(fig.add_subplot(2, 2, 1, projection="3d"),
               logQ2, logK2, g1, sigma_a, "g1", r"$\Sigma_A$", 1)
    _plot_pair(fig.add_subplot(2, 2, 2, projection="3d"),
               logQ2, logK2, g2, 2 * delta_a, "g2", r"$2\Delta_A$", 2)
    _plot_pair(fig.add_subplot(2, 2, 3, projection="3d"),
               logQ2, logK2, g3, -2 * delta_b, "g3", r"$-2\Delta_B$", 3)
    _plot_pair(fig.add_subplot(2, 2, 4, projection="3d"),
               logQ2, logK2, g4, None, "g4", "", 4)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(outdir, "wti.pdf")
    fig.savefig(out)
    plt.close(fig)
    return out
