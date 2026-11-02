"""qpv — reader and plotting utilities for the Quark_Photon_Vertex HDF5 output.

Replaces the old gnuplot_script/*.gplt + plot_hvp.py pipeline. See reader.py
for loading a run and plot_*.py for the individual figures.
"""

from .reader import Run, open_run

__all__ = ["Run", "open_run"]
