#!/usr/bin/env bash
# Run the quark-photon-vertex solve for several flavours, each into its own
# HDF5 file. Non-light flavours require the quark DSE (-d), so 'd' is always in
# the flag set here.
#
# Each production run takes ~1-2 h; this script runs them sequentially in the
# foreground. To parallelise, launch individual invocations in the background
# instead. Outputs land in diag_runs/<flavour>_<mode>/output_<...>.h5.
#
# Usage:
#   ./run_flavors.sh [BINARY] [FLAGS] [FLAVOURS...]
# Examples:
#   ./run_flavors.sh                                  # all four, flags "dphv"
#   ./run_flavors.sh build/Quark_Photon_Vertex/Quark_Photon_Vertex dp l s
set -euo pipefail

BIN="${1:-build/Quark_Photon_Vertex/Quark_Photon_Vertex}"
FLAGS="${2:-dphv}" # d = DSE (required for s/c/b), p = Pauli-Villars, h = HVP
shift || true
shift || true
FLAVOURS=("${@:-l s c b}")
# handle the default (unquoted expansion above yields one element "l s c b")
if [[ ${#FLAVOURS[@]} -eq 1 && "${FLAVOURS[0]}" == "l s c b" ]]; then
  FLAVOURS=(l s c b)
fi

if [[ ! -x "$BIN" ]]; then
  echo "error: binary not found or not executable: $BIN" >&2
  exit 1
fi
BIN="$(readlink -f "$BIN")"

for fl in "${FLAVOURS[@]}"; do
  tag="${fl}_${FLAGS}"
  outdir="diag_runs/${tag}"
  mkdir -p "$outdir"
  echo "=== flavour '${fl}' (flags '${FLAGS}${fl}') -> ${outdir}/ ==="
  (cd "$outdir" && "$BIN" "${FLAGS}${fl}")
done

echo "Done. HDF5 files are under diag_runs/<flavour>_<mode>/."
