#!/usr/bin/env bash
# Tolerance-aware diff of bench results against the frozen baseline.
#
# Usage:
#   bench/diff_baseline.sh                     # diff bench/results/current
#   bench/diff_baseline.sh path/to/results     # diff arbitrary results dir
#   RTOL=1e-6 ATOL=1e-10 bench/diff_baseline.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS_DIR="${1:-${REPO_ROOT}/bench/results/current}"
BASELINE_DIR="${REPO_ROOT}/tests/bench_baseline"
RTOL="${RTOL:-1e-6}"
# atol = 1e-7 is 100× tighter than the iteration's own target_acc (1e-5),
# strict enough to catch real regressions, lenient enough to absorb the FP
# rounding noise that build flags like -march=native or LTO introduce.
ATOL="${ATOL:-1e-7}"
COMBOS_DEFAULT="nodse_nopv nodse_pv dse_nopv dse_pv"
COMBOS="${COMBOS:-$COMBOS_DEFAULT}"

# Files excluded from regression diffing. Delta_A = (A(k+²) - A(k-²)) /
# (k+² - k-²) suffers intrinsic catastrophic cancellation when q² is small
# and z is near ±1; any FP-affecting change (build flags, OpenMP schedule,
# DSE convergence path) drifts these values by O(1e-3) relative. The WTI
# is post-processing only — it does not feed back into the BSE iteration.
EXCLUDE_PATTERN='^(w_file_idx_1|w_z0_file_idx_1)\.dat$'

if [ ! -d "${BASELINE_DIR}" ] || \
   [ -z "$(find "${BASELINE_DIR}" -mindepth 2 -name '*.dat' -print -quit 2>/dev/null)" ]; then
    echo "No baseline at ${BASELINE_DIR}." >&2
    echo "Capture one first:  bench/run_bench.sh ${BASELINE_DIR}" >&2
    exit 2
fi

failed=0
checked=0
for combo in ${COMBOS}; do
    combo_dir="${BASELINE_DIR}/${combo}/"
    rdir="${RESULTS_DIR}/${combo}"
    if [ ! -d "${combo_dir}" ]; then
        echo "skipping ${combo}: no baseline subdir"
        continue
    fi
    if [ ! -d "${rdir}" ]; then
        echo "missing combo dir: ${rdir}"
        failed=1
        continue
    fi
    for f in "${combo_dir}"*.dat; do
        [ -f "$f" ] || continue
        name="$(basename "$f")"
        if [[ "${name}" =~ ${EXCLUDE_PATTERN} ]]; then
            echo "skip  ${combo}/${name} (WTI catastrophic cancellation)"
            continue
        fi
        rfile="${rdir}/${name}"
        if [ ! -f "${rfile}" ]; then
            echo "FAIL ${combo}/${name}: missing in results"
            failed=1
            continue
        fi
        if ! python3 "${REPO_ROOT}/bench/numdiff.py" \
                --rtol "${RTOL}" --atol "${ATOL}" "$f" "${rfile}"; then
            failed=1
        fi
        checked=$((checked + 1))
    done
done

echo "---"
if [ $failed -ne 0 ]; then
    echo "DIFF FAILED (${checked} files compared, rtol=${RTOL} atol=${ATOL})"
    exit 1
fi
echo "OK: ${checked} files match within rtol=${RTOL} atol=${ATOL}"
