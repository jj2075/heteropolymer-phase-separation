#!/bin/bash
set -euf -o pipefail

# Detect GPU count: fallback to 1 if no GPUs detected
gpu_count=${1:-$(nvidia-smi --list-gpus | wc -l || echo 1)}

echo "Running simulation on ${gpu_count} GPU(s)..."

# Run LAMMPS with Kokkos (GPU-accelerated)
mpirun -n ${gpu_count} lmp -k on g ${gpu_count} -sf kk -pk kokkos cuda/aware on neigh full comm device binsize 2.8 -in equilibrate.in
