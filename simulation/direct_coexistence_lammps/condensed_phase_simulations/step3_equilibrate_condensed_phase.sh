#!/bin/bash
#SBATCH --job-name=equilib       # job name for longer NVT equilibration
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks=1               # single task
#SBATCH --cpus-per-task=1        # number of CPU cores per task
#SBATCH --mem-per-cpu=800M       # memory per CPU core
#SBATCH --gres=gpu:1             # number of GPUs per node
#SBATCH --time=18:59:00          # time limit (HH:MM:SS)
#SBATCH --array=1-30             # array of jobs for multiple sequences

set -euf -o pipefail

# This script runs a longer NVT equilibration step for each sequence after temperature cycling

seq_index=${SLURM_ARRAY_TASK_ID}
data_dir="<path-to-data-files>"  # dir containing epsilon and sequence data
epsilon_file="$data_dir/epsilon.txt"  # file with seqID and epsilon values
seqs_file="$data_dir/seqs-ps.lst"  # file with sequence IDs and sequences

# Get seqID, sequence, and epsilon for the current job
sequence_line=$(sed -n "${seq_index}p" "$seqs_file")
seqid=$(echo "$sequence_line" | awk '{print $1}')
sequence=$(echo "$sequence_line" | awk '{print $2}')
epsilon=$(awk -v id="$seqid" '$1 == id {print $2}' "$epsilon_file")

echo "Sequence ID: $seqid"
echo "Sequence: $sequence"
echo "Epsilon: $epsilon"

# Create dir for this sequence
mkdir -p "$seqid"
cd "$seqid"

# Copy input files for longer NVT equilibration
keyfiles_dir="<path-to-keyfiles>"  # dir with keyfiles for NVT equilibration
cp "$keyfiles_dir/run_equilibrate.sh" .
cp "$keyfiles_dir/equilibrate.in" .

# Replace placeholders in input files with sequence-specific values
sed -i "s/{epsilon}/$epsilon/g" equilibrate.in
sed -i "s/{trial_number}/1/g" equilibrate.in

module purge
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
 
singularity run --nv -B $PWD:/host_pwd --pwd /host_pwd $HOME/software/lammps_container/lammps_10Feb2021.sif ./run_equilibrate.sh > "${seqid}-equilibrate-trial1.log" 2>&1

