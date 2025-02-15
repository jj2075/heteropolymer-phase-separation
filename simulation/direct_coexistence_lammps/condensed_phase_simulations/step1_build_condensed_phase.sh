#!/bin/bash
#SBATCH --job-name=cond_phase_init       # Job name
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1                   # Total number of tasks
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem-per-cpu=800M           # Memory allocation per CPU core
#SBATCH --gres=gpu:1                 # Request 1 GPU
#SBATCH --time=02:30:00              # Time limit (HH:MM:SS)
#SBATCH --array=1-72                 # Array of jobs for 72 sequences

# Generalized and reproducible job script for direct coexistence simulations.
# This script initializes and runs phase compression for each sequence.

set -euf -o pipefail

line_number=${SLURM_ARRAY_TASK_ID}
data_dir="<path-to-data-files>"  # dir containing e psilon and sequence data
epsilon_file="$data_dir/epsilon.csv"  # file with seqID and corresponding epsilon values
seqs_file="$data_dir/seqs-ps.lst"  # file with sequence IDs and sequences

# Get seqID, sequence, and epsilon for the current job from CSV file
subdir_line=$(sed -n "${line_number}p" "$epsilon_file")
seqid=$(echo "$subdir_line" | awk -F ',' '{print $1}')      # seqid in the first column
sequence=$(echo "$subdir_line" | awk -F ',' '{print $2}')  # sequence in the second column
epsilon=$(echo "$subdir_line" | awk -F ',' '{print $4}')   # epsilon value in the fourth column

echo "Sequence ID: $seqid"
echo "Sequence: $sequence"
echo "Epsilon: $epsilon"

# Create dir for this sequence
mkdir -p "$seqid"
cd "$seqid"

# Copy input templates for phase compression
keyfiles_dir="<path-to-keyfiles>"  # Directory with keyfiles for initialization
cp "$keyfiles_dir/compress-0.in" .
cp "$keyfiles_dir/phase.json" .
cp "$keyfiles_dir/run_compress.sh" .

# replace epsilon and sequence values
sed -i "s/{epsilon}/$epsilon/g" compress-0.in
sed -i "s/{sequence}/$sequence/g" phase.json

# Generate the initial configuration
python3 "$keyfiles_dir/gen_init_config.py" phase.json phase-

module purge
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Run the simulation using sbatch
# Example:
# srun ./run_compress.sh
