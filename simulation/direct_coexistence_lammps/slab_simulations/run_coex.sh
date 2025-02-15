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

seq_index=${SLURM_ARRAY_TASK_ID}
data_dir="<path-to-data-files>"  # dir containing epsilon and sequence data
epsilon_file="$data_dir/epsilon.csv"  # file with seqID and epsilon values
seqs_file="$data_dir/seqs-ps.lst"  # file with sequence IDs and sequences

# Get seqID, sequence, and epsilon for the current job from CSV file
subdir_line=$(sed -n "${line_number}p" "$epsilon_file")
seqid=$(echo "$subdir_line" | awk -F ',' '{print $1}')      # seqid in the first column
sequence=$(echo "$subdir_line" | awk -F ',' '{print $2}')  # sequence in the second column
epsilon=$(echo "$subdir_line" | awk -F ',' '{print $4}')   # epsilon value in the fourth column

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
