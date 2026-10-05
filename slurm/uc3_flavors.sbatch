#!/bin/bash
# All four quark flavours (light, strange, charm, bottom) with the production
# grid from parameters.hh, on bwUniCluster 3.0. One array task per flavour.
#
# Submit from the repository root:
#   sbatch slurm/uc3_flavors.sbatch
# Outputs: results/<flavour>/output_<flavour>_dse_pv.h5 (+ run.log)
#
# Each task builds its own binary in node-local $TMPDIR: CMakeLists.txt uses
# -march=native, so the build must happen on the node type that runs it, and
# separate build dirs avoid races between array tasks.
#
#SBATCH --job-name=qpv_flavors
#SBATCH --partition=cpu
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --threads-per-core=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=slurm-qpv-%A_%a.out

set -euo pipefail

FLAVOURS=(l s c b)
NAMES=(light strange charm bottom)
fl=${FLAVOURS[$SLURM_ARRAY_TASK_ID]}
name=${NAMES[$SLURM_ARRAY_TASK_ID]}

module purge
module load compiler/gnu/14.2
module load lib/hdf5/1.14-gnu-14.2-serial

src=$SLURM_SUBMIT_DIR
build=${TMPDIR:-/tmp}/qpv_build_${SLURM_JOB_ID}
outdir=$src/results/$name
mkdir -p "$outdir"

cmake -S "$src" -B "$build" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=g++ -DBUILD_TESTS=OFF > "$outdir/build.log" 2>&1
cmake --build "$build" -j "$SLURM_CPUS_PER_TASK" --target Quark_Photon_Vertex \
      >> "$outdir/build.log" 2>&1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "flavour $name on $(hostname), $OMP_NUM_THREADS threads, $(date)"
cd "$outdir"
# d = quark DSE (required for s/c/b), p = Pauli-Villars, h = HVP
srun "$build/Quark_Photon_Vertex/Quark_Photon_Vertex" "dph$fl" > run.log 2>&1
echo "done $(date)"
