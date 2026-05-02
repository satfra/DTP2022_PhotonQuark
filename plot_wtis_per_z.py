#!/usr/bin/env python3
"""Plot WTI residuals at every z-grid point (no z=0 averaging).

Reads the full 3D output: w_file_idx_{0,1,2}.dat (Sigma_A, Delta_A, Delta_B) and
fg_file_idx_{8,9,10,11}.dat (g1..g4). Produces, per WTI:

  1. A heatmap of |g_i - target_i| over (z, k^2) at a representative Q^2.
  2. A line plot of max-over-k relative residual vs z, one curve per Q^2.

The goal is to check whether the WTI residual is concentrated at the innermost
|z| nodes (consistent with the 1/(z*s) basis-amplification diagnosis) or spread
uniformly across z (which would point at a different cause).

Usage: python plot_wtis_per_z.py <folder> [--out NAME.pdf]
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_3d(path: Path):
    """Return (q_axis, k_axis, z_axis, V[q,k,z]) from a w_file_idx_* / fg_file_idx_* file."""
    data = np.loadtxt(path, comments="#")
    q, _i, k, z, re = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    q_axis = np.unique(q)
    k_axis = np.unique(k)
    z_axis = np.unique(z)
    nq, nk, nz = len(q_axis), len(k_axis), len(z_axis)
    if nq * nk * nz != len(re):
        raise ValueError(f"{path}: {len(re)} rows do not factor as {nq}*{nk}*{nz}")
    # Saved as q outer → k → z inner (see fileIO.hh ordering).
    V = re.reshape(nq, nk, nz)
    return q_axis, k_axis, z_axis, V


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", type=Path)
    ap.add_argument("--out", default="wti_per_z.pdf")
    args = ap.parse_args()

    folder: Path = args.folder
    needed = [folder / f"fg_file_idx_{i}.dat" for i in (8, 9, 10, 11)] + \
             [folder / f"w_file_idx_{i}.dat" for i in (0, 1, 2)]
    missing = [p for p in needed if not p.is_file()]
    if missing:
        print(f"error: missing files in {folder}:", file=sys.stderr)
        for p in missing:
            print(f"  {p.name}", file=sys.stderr)
        sys.exit(1)

    q_axis, k_axis, z_axis, g1 = load_3d(folder / "fg_file_idx_8.dat")
    _, _, _, g2 = load_3d(folder / "fg_file_idx_9.dat")
    _, _, _, g3 = load_3d(folder / "fg_file_idx_10.dat")
    _, _, _, g4 = load_3d(folder / "fg_file_idx_11.dat")
    _, _, _, sig_A = load_3d(folder / "w_file_idx_0.dat")
    _, _, _, del_A = load_3d(folder / "w_file_idx_1.dat")
    _, _, _, del_B = load_3d(folder / "w_file_idx_2.dat")

    # Per-WTI residual arrays, full 3D (q, k, z).
    res = {
        1: ("g1 vs Sigma_A",       g1, sig_A),
        2: ("g2 vs 2*Delta_A",     g2, 2.0 * del_A),
        3: ("g3 vs -2*Delta_B",    g3, -2.0 * del_B),
        4: ("g4 (target = 0)",     g4, np.zeros_like(g4)),
    }

    fig, axes = plt.subplots(4, 2, figsize=(13, 16))
    fig.suptitle(f"WTI per-z diagnostic — {folder.name}", fontsize=13)

    for row, (n, (label, X, T)) in enumerate(res.items()):
        diff = X - T
        peak = np.maximum(np.abs(T).max(axis=(1,)), np.abs(X).max(axis=(1,)) * 1e-12)
        # max-over-k absolute residual per (q, z), normalised by peak |T|+|X|.
        peak_qz = np.maximum(np.abs(T).max(axis=1, keepdims=False),
                             np.abs(X).max(axis=1, keepdims=False) * 1e-12)
        rel_qz = np.abs(diff).max(axis=1) / np.maximum(peak_qz, 1e-300)

        # --- Left: line plot, one curve per q, residual vs z. ---
        ax = axes[row, 0]
        for qi, q_val in enumerate(q_axis):
            ax.plot(z_axis, rel_qz[qi, :], lw=0.8, alpha=0.6,
                    label=f"Q²={q_val:.1e}" if qi in (0, len(q_axis) // 2, len(q_axis) - 1) else None)
        ax.set_xlabel("z")
        ax.set_ylabel("max_k |g - target| / max |target|")
        ax.set_yscale("log")
        ax.set_title(f"WTI {n}: {label}\n(curves = different Q²)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper center")

        # --- Right: heatmap of relative residual over (z, q), at the worst-k slice ---
        ax = axes[row, 1]
        # rel residual at the worst-k for each (q,z): use rel_qz directly.
        Z, Q2 = np.meshgrid(z_axis, q_axis)
        # Use log10 of relative residual (clip at 1e-12 to avoid -inf).
        log_rel = np.log10(np.maximum(rel_qz, 1e-12))
        im = ax.pcolormesh(Z, Q2, log_rel, shading="auto", cmap="viridis")
        ax.set_yscale("log")
        ax.set_xlabel("z")
        ax.set_ylabel("Q² [GeV²]")
        ax.set_title(f"WTI {n}: log₁₀ max_k rel.err.")
        plt.colorbar(im, ax=ax)

        global_max = float(np.nanmax(rel_qz))
        ax.text(0.02, 0.98, f"global max = {global_max:.2e}",
                transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = folder / args.out
    fig.savefig(out)
    print(f"wrote {out}")
    plt.close(fig)

    # Print a concise summary line per WTI.
    print("\nSummary (max-over-(q,k,z) of relative residual):")
    for n, (label, X, T) in res.items():
        diff = X - T
        peak = max(float(np.abs(T).max()), 1e-300)
        m_abs = float(np.abs(diff).max())
        rel = m_abs / peak
        print(f"  WTI{n}: {label:<28s}  max|abs|={m_abs:.3e}   max rel={rel:.3e}")


if __name__ == "__main__":
    main()
