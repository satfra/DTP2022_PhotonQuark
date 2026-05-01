#!/usr/bin/env bash
# Build with -DBENCHMARK_PROFILE=ON and run all four flag combos, capturing
# .dat outputs and per-combo wall time into the chosen output directory.
#
# Usage:
#   bench/run_bench.sh                                 # → bench/results/current
#   bench/run_bench.sh tests/bench_baseline             # capture baseline
#   COMBOS="nodse_nopv dse_pv" bench/run_bench.sh       # subset
#
# Optionally tag the build with TIER=tier1 etc; this only affects logging.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_BASE="${1:-${REPO_ROOT}/bench/results/current}"
BUILD_DIR="${REPO_ROOT}/build_bench"
TIER_TAG="${TIER:-current}"
COMBOS_DEFAULT="nodse_nopv nodse_pv dse_nopv dse_pv"
COMBOS="${COMBOS:-$COMBOS_DEFAULT}"

# Map result-dir name → flag string for the executable.
flags_for() {
    case "$1" in
        nodse_nopv) echo "" ;;
        nodse_pv)   echo "p" ;;
        dse_nopv)   echo "d" ;;
        dse_pv)     echo "dp" ;;
        *) echo "unknown combo: $1" >&2; exit 2 ;;
    esac
}

mkdir -p "${OUT_BASE}"
RESULTS_LOG="${OUT_BASE}/timings.txt"
: > "${RESULTS_LOG}"

echo ">>> Configuring (BENCHMARK_PROFILE=ON, Release)..."
cmake -B "${BUILD_DIR}" -S "${REPO_ROOT}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBENCHMARK_PROFILE=ON > "${OUT_BASE}/cmake.log" 2>&1
echo ">>> Building Quark_Photon_Vertex..."
cmake --build "${BUILD_DIR}" -j --target Quark_Photon_Vertex \
    > "${OUT_BASE}/build.log" 2>&1

EXE="${BUILD_DIR}/Quark_Photon_Vertex/Quark_Photon_Vertex"
[ -x "${EXE}" ] || { echo "missing executable: ${EXE}" >&2; exit 1; }

echo "tier=${TIER_TAG}" >> "${RESULTS_LOG}"
echo "host=$(hostname) cores=$(nproc) date=$(date -Iseconds)" >> "${RESULTS_LOG}"

for combo in ${COMBOS}; do
    flags="$(flags_for "${combo}")"
    out="${OUT_BASE}/${combo}"
    rm -rf "${out}"
    mkdir -p "${out}"
    echo ">>> Running ${combo} (flags='${flags}')..."

    pushd "${out}" > /dev/null
    start_ns=$(date +%s%N)
    if [ -z "${flags}" ]; then
        "${EXE}" > "${combo}.stdout" 2> "${combo}.stderr"
    else
        "${EXE}" "${flags}" > "${combo}.stdout" 2> "${combo}.stderr"
    fi
    end_ns=$(date +%s%N)
    popd > /dev/null

    elapsed_ms=$(( (end_ns - start_ns) / 1000000 ))
    elapsed_s=$(awk -v ms="${elapsed_ms}" 'BEGIN{printf "%.3f", ms/1000.}')
    echo "${combo}: ${elapsed_s}s" | tee -a "${RESULTS_LOG}"
done

echo ">>> Done. Outputs under ${OUT_BASE}/, timings in ${RESULTS_LOG}"
