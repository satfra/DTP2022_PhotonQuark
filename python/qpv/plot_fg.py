"""Surface plots of the z-averaged f/g dressing functions.

Reproduces gnuplot_script/fg.gplt: for each of f1..f8, g1..g4 a 3D surface of
Re(fg_z0) over (log10 Q^2, log10 k^2). The old gnuplot used `splot ... u 1:3:4`
on fg_z0_file, i.e. q_sq : k_sq : Re(fg).
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .reader import Run  # noqa: E402

# f1,f2,f5..f8 were z-clamped to [-2,2] in the old gnuplot; keep that for
# comparability, leave f3,f4,g1..g4 auto-scaled.
_ZCLAMP = {"f1", "f2", "f5", "f6", "f7", "f8"}


def _surface(ax, run: Run, name: str) -> None:
    data = np.real(run.fgz0(name))          # [q, k]
    logq = np.log10(run.q_grid)
    logk = np.log10(run.k_sq)
    K, Q = np.meshgrid(logk, logq)          # both [q, k]
    ax.plot_surface(Q, K, data, cmap="viridis", linewidth=0, antialiased=True)
    ax.set_xlabel(r"$\log_{10} Q^2$")
    ax.set_ylabel(r"$\log_{10} k^2$")
    ax.set_zlabel(rf"$\mathrm{{Re}}\,{name}$")
    ax.set_title(f"{name}  —  {run.label}")
    if name in _ZCLAMP:
        ax.set_zlim(-2, 2)


def plot_one(run: Run, name: str, outdir: str) -> str:
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    _surface(ax, run, name)
    path = os.path.join(outdir, f"fg_{name}.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_all_fg(run: Run, outdir: str) -> list[str]:
    os.makedirs(outdir, exist_ok=True)
    return [plot_one(run, name, outdir) for name in run.fg_names]
