#!/bin/bash
#SBATCH --job-name=direct-coex    # Job name
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Single task
#SBATCH --cpus-per-task=1         # CPU cores per task
#SBATCH --mem-per-cpu=840M        # Memory per CPU core
#SBATCH --gres=gpu:1              # Number of GPUs per node
#SBATCH --time=23:59:00           # Total run time limit (HH:MM:SS)
#SBATCH --array=26-27             # Array range for sequences

set -euf -o pipefail

data_dir="<path-to-data-files>"  # dir containing sequence and epsilon data
epsilon_file="$data_dir/epsilon.txt"
seqs_file="$data_dir/seqs-ps.lst"

# Get seqID and epsilon for the current job
line_number=${SLURM_ARRAY_TASK_ID}
subdir_line=$(sed -n "${line_number}p" "$epsilon_file")
seqid=$(echo "$subdir_line" | awk '{print $1}')
epsilon=$(awk -v id="$seqid" '$1 == id {print $2}' "$epsilon_file")
sequence=$(awk -v id="$seqid" '$1 == id {print $2}' "$seqs_file")

echo "Sequence ID: $seqid"
echo "Sequence: $sequence"
echo "Epsilon: $epsilon"

# setup dir for the current sequence
mkdir -p "$seqid"
keyfiles_dir="<path-to-keyfiles>"
cp $keyfiles_dir/coex.in "$seqid/"
cp $keyfiles_dir/run_lammps.sh "$seqid/"

# Replace with seq-specific B2-matched epsilon value
sed -i "s/{epsilon}/$epsilon/g" "$seqid/coex.in"

module purge
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$seqid"
singularity run --nv -B $PWD:/host_pwd --pwd /host_pwd $HOME/software/lammps_container/lammps_10Feb2021.sif ./run_lammps.sh > "${seqid}.log" 2>&1
wait

echo "Computing density profile..."
# Analyze density profiles using the last 1000 snapshots
# coex_zbincounts.dat contains density profile data every 20,000 timesteps
python3 $keyfiles_dir/analyze_zbins.py coex_zbincounts.dat --neq 23000
python3 $keyfiles_dir/plot_density_profile.py --seqid "$seqid" --epsilon "$epsilon"
