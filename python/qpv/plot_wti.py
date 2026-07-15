"""Ward-Takahashi-identity comparison plots.

Reproduces gnuplot_script/wti.gplt. In the exact-WTI limit the g-functions
match the propagator quantities:
    WTI1:  g1  ==  Sigma_A
    WTI2: -g2  ==  2 * Delta_A
    WTI3:  g3  == -2 * Delta_B
    WTI4:  g4  ==  0
Each panel overlays the two z-averaged surfaces over (log10 Q^2, log10 k^2).
(The physically correct Delta_B index is used for WTI3; the original gnuplot
mislabelled it via w_z0 idx 0.)
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .reader import Run  # noqa: E402

# (title, g-name, g-scale, reference-callable-or-None, reference-label)
def _specs(run: Run):
    return [
        ("wti1", "g1", 1.0, lambda: np.real(run.wz0("Sigma_A")),      r"$\Sigma_A$"),
        ("wti2", "g2", -1.0, lambda: 2.0 * np.real(run.wz0("Delta_A")), r"$2\Delta_A$"),
        ("wti3", "g3", 1.0, lambda: -2.0 * np.real(run.wz0("Delta_B")), r"$-2\Delta_B$"),
        ("wti4", "g4", 1.0, None, None),
    ]


def _panel(run: Run, tag: str, gname: str, gscale: float, ref, ref_label: str, outdir: str) -> str:
    logq = np.log10(run.q_grid)
    logk = np.log10(run.k_sq)
    K, Q = np.meshgrid(logk, logq)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")

    gdata = gscale * np.real(run.fgz0(gname))
    glabel = (r"$-%s$" % gname) if gscale < 0 else r"$%s$" % gname
    ax.plot_surface(Q, K, gdata, color="tab:blue", alpha=0.7, linewidth=0, label=glabel)

    handles = [plt.Line2D([0], [0], color="tab:blue", label=glabel)]
    if ref is not None:
        ax.plot_surface(Q, K, ref(), color="tab:orange", alpha=0.5, linewidth=0)
        handles.append(plt.Line2D([0], [0], color="tab:orange", label=ref_label))

    ax.set_xlabel(r"$\log_{10} Q^2$")
    ax.set_ylabel(r"$\log_{10} k^2$")
    ax.set_zlabel("value")
    ax.set_title(f"{tag}: {glabel}" + (f" vs {ref_label}" if ref is not None else "")
                 + f"  —  {run.label}")
    ax.legend(handles=handles, loc="upper right")

    path = os.path.join(outdir, f"{tag}.pdf")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_all_wti(run: Run, outdir: str) -> list[str]:
    os.makedirs(outdir, exist_ok=True)
    return [_panel(run, tag, g, s, ref, lbl, outdir) for tag, g, s, ref, lbl in _specs(run)]
