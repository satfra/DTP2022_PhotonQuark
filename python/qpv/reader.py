"""Load a Quark_Photon_Vertex HDF5 run into numpy arrays.

The C++ side (hdf5IO.hh) writes complex data as a compound {r, i} datatype,
which h5py maps to numpy complex128 automatically — so fg/b/w come back as
native complex arrays with q as the leading axis.
"""

from __future__ import annotations

from dataclasses import dataclass

import h5py
import numpy as np


def _decode_names(ds) -> list[str]:
    return [s.decode() if isinstance(s, bytes) else str(s) for s in ds[:]]


@dataclass
class Run:
    """One HDF5 run. Complex tensors have q as the leading axis."""

    path: str
    # grids
    q_grid: np.ndarray      # [q]  external Q^2
    k_sq: np.ndarray        # [k]  physical k^2
    k_grid: np.ndarray      # [k]  raw log(k^2)
    z_grid: np.ndarray      # [z]
    y_grid: np.ndarray      # [y]
    # result tensors
    fg: np.ndarray          # [q,12,k,z] complex   f1..f8, g1..g4
    fg_z0: np.ndarray       # [q,12,k]   complex   z-averaged f/g
    b: np.ndarray           # [q,12,k,z] complex
    w: np.ndarray           # [q,3,k,z]  complex   Sigma_A, Delta_A, Delta_B
    w_z0: np.ndarray        # [q,3,k]    complex
    # optional HVP
    hvp_p_sq: np.ndarray | None
    hvp_pi: np.ndarray | None
    # optional quark-DSE solution (present only for -d runs)
    dse: dict | None       # keys: p_sq, log_p_sq, A, B, M (each [n_dse])
    # metadata
    attrs: dict
    fg_names: list[str]
    w_names: list[str]

    # --- named-slice helpers -------------------------------------------------
    def fg_index(self, name: str) -> int:
        return self.fg_names.index(name)

    def w_index(self, name: str) -> int:
        return self.w_names.index(name)

    def fgz0(self, name: str) -> np.ndarray:
        """z-averaged f/g structure by name, shape [q, k]."""
        return self.fg_z0[:, self.fg_index(name), :]

    def wz0(self, name: str) -> np.ndarray:
        """z-averaged WTI quantity by name, shape [q, k]."""
        return self.w_z0[:, self.w_index(name), :]

    @property
    def has_hvp(self) -> bool:
        return self.hvp_p_sq is not None

    @property
    def has_dse(self) -> bool:
        return self.dse is not None

    @property
    def label(self) -> str:
        """Human-readable run tag from metadata (flavour + mode)."""
        tag = str(self.attrs.get("flavor", "?"))
        if self.attrs.get("use_dse"):
            tag += " (DSE"
            tag += ", PV)" if self.attrs.get("use_pauli_villars") else ")"
        elif self.attrs.get("use_pauli_villars"):
            tag += " (PV)"
        return tag


def _legacy_pi_hat(p_sq: np.ndarray, pi_old: np.ndarray, z2: float,
                   fit_max_p_sq: float = 1e-2) -> np.ndarray:
    """Pi(p^2) - Pi(0) from a file written before hvp.hh stored it directly.

    Those files hold half the traced loop T (the y-angle factor 2 was missing
    from the measure), subtracted at p_sq[0]. Mirrors hvp::hvp_driver: fit
    T = C + s p^2 + O(p^4) at low p^2, then Pi_hat = z2/3 (s - (T - C)/p^2).
    """
    trace = 2.0 * pi_old
    n = max(int(np.sum(p_sq <= fit_max_p_sq)), 6)
    C, s = np.polynomial.polynomial.polyfit(p_sq[:n], trace[:n], 3)[:2]
    return z2 / 3.0 * (s - (trace - C) / p_sq)


def open_run(path: str) -> Run:
    with h5py.File(path, "r") as f:
        g = f["grids"]
        hvp_p = f["hvp/p_sq"][:] if "hvp" in f else None
        hvp_pi = f["hvp/Pi"][:] if "hvp" in f else None
        if hvp_pi is not None and "pi_at_zero" not in f["hvp"].attrs:
            hvp_pi = _legacy_pi_hat(hvp_p, hvp_pi, float(f.attrs.get("z2", 1.0)))
        dse = None
        if "quark_dse" in f:
            gd = f["quark_dse"]
            dse = {k: gd[k][:] for k in ("p_sq", "log_p_sq", "A", "B", "M")}
        return Run(
            path=path,
            q_grid=g["q_grid"][:],
            k_sq=g["k_sq"][:],
            k_grid=g["k_grid"][:],
            z_grid=g["z_grid"][:],
            y_grid=g["y_grid"][:],
            fg=f["fg"][:],
            fg_z0=f["fg_z0"][:],
            b=f["b"][:],
            w=f["w"][:],
            w_z0=f["w_z0"][:],
            hvp_p_sq=hvp_p,
            hvp_pi=hvp_pi,
            dse=dse,
            attrs={k: v for k, v in f.attrs.items()},
            fg_names=_decode_names(f["fg_index_names"]),
            w_names=_decode_names(f["w_index_names"]),
        )
