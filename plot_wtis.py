#!/usr/bin/env python3
"""Plot the four Ward-Takahashi identity checks for a quark-photon vertex run.

Usage: python plot_wtis.py <folder>

The folder must contain the seven z=0 output files emitted by the simulation:
  fg_z0_file_idx_{8,9,10,11}.dat   (g1..g4)
  w_z0_file_idx_{0,1,2}.dat        (Sigma_A, Delta_A, Delta_B)
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)


def load_grid(path: Path):
    """Return (Q2, K2, Z) meshgrids reshaped from the q_sq, i, k_sq, Re, Im rows."""
    data = np.loadtxt(path, comments="#")
    q2, k2, re = data[:, 0], data[:, 2], data[:, 3]
    q_axis = np.unique(q2)
    k_axis = np.unique(k2)
    n_q, n_k = len(q_axis), len(k_axis)
    if n_q * n_k != len(re):
        raise ValueError(f"{path}: {len(re)} rows do not factor as {n_q}*{n_k}")
    Z = re.reshape(n_q, n_k)
    Q2, K2 = np.meshgrid(q_axis, k_axis, indexing="ij")
    return Q2, K2, Z


def plot_pair(ax, Q2, K2, Zg, Zt, g_label, t_label, wti_num):
    logQ2, logK2 = np.log10(Q2), np.log10(K2)
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
        # L_inf relative error: max|g - target| / max|target|. Pointwise |g-t|/|t|
        # blows up wherever target is near a zero crossing while g is small but
        # nonzero — those points are visually invisible but numerically inflate
        # the error. Peak-normalising matches what the overlaid surfaces show.
        max_diff = float(np.abs(Zg - Zt).max())
        peak = float(np.abs(Zt).max())
        max_rel = max_diff / peak if peak > 0 else float("nan")
        info = f"max rel. err. = {max_rel:.2e}"
    else:
        info = f"max |{g_label}| = {float(np.nanmax(np.abs(Zg))):.2e}"
    ax.text2D(0.95, 0.95, info, transform=ax.transAxes, fontsize=9,
              ha="right", va="top",
              bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, help="folder containing the .dat files")
    args = parser.parse_args()

    folder: Path = args.folder
    needed = [folder / f"fg_z0_file_idx_{i}.dat" for i in (8, 9, 10, 11)] + \
             [folder / f"w_z0_file_idx_{i}.dat" for i in (0, 1, 2)]
    missing = [p for p in needed if not p.is_file()]
    if missing:
        print(f"error: missing files in {folder}:", file=sys.stderr)
        for p in missing:
            print(f"  {p.name}", file=sys.stderr)
        sys.exit(1)

    Q2, K2, g1 = load_grid(folder / "fg_z0_file_idx_8.dat")
    _, _, g2 = load_grid(folder / "fg_z0_file_idx_9.dat")
    _, _, g3 = load_grid(folder / "fg_z0_file_idx_10.dat")
    _, _, g4 = load_grid(folder / "fg_z0_file_idx_11.dat")
    _, _, sigma_a = load_grid(folder / "w_z0_file_idx_0.dat")
    _, _, delta_a = load_grid(folder / "w_z0_file_idx_1.dat")
    _, _, delta_b = load_grid(folder / "w_z0_file_idx_2.dat")

    fig = plt.figure(figsize=(14, 11))
    fig.suptitle(f"WTI checks — {folder.name}", fontsize=14)

    plot_pair(fig.add_subplot(2, 2, 1, projection="3d"),
              Q2, K2, g1, sigma_a, "g1", r"$\Sigma_A$", 1)
    plot_pair(fig.add_subplot(2, 2, 2, projection="3d"),
              Q2, K2, g2, 2 * delta_a, "g2", r"$2\Delta_A$", 2)
    plot_pair(fig.add_subplot(2, 2, 3, projection="3d"),
              Q2, K2, g3, -2 * delta_b, "g3", r"$-2\Delta_B$", 3)
    plot_pair(fig.add_subplot(2, 2, 4, projection="3d"),
              Q2, K2, g4, None, "g4", "", 4)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = folder / "wti.pdf"
    fig.savefig(out)
    print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
