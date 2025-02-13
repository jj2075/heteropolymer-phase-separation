#!/bin/bash
#SBATCH --job-name=rg_sim          # Job name
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --cpus-per-task=1         # Number of CPU cores per task
#SBATCH --mem-per-cpu=250M        # Memory allocation per CPU core
#SBATCH --time=01:30:00           # Time limit (HH:MM:SS)
#SBATCH --array=1-72              # Array of 72 jobs (each for a different sequence)
#SBATCH --constraint=cascade,skylake  # Architecture constraints

# This job script runs LAMMPS simulations to compute Rg over time for a single B2-matched heteropolymer chain.
# The epsilon value is derived from B2-matching and passed into the simulation.

ulimit -s unlimited
ulimit -c unlimited

module purge
module load intel/2021.1
module load intel-mpi/gcc/2021.1.1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

line_number=$((SLURM_ARRAY_TASK_ID))

# Generalized paths
data_dir="<path_to_data_files>" # Dir containing data files for sequences and epsilon values
epsilon_file="$data_dir/combined_epsilon_1000.txt"  # Contains seqID and corresponding epsilon values
seqs_file="$data_dir/b2_1000_combined_seqs-ps.lst"  # Contains sequence IDs
current_simulation_rg_dir="<path-to-simulation_scripts>/rg"  # Directory for Rg simulations

line=$(sed -n "${line_number}p" "$epsilon_file")
seqid=$(echo "$line" | awk '{print $1}')
epsilon=$(echo "$line" | awk '{print $2}')
sequence=$(awk -v id="$seqid" '$1 == id {print $2}' "$seqs_file")

echo "seqID: $seqid"
echo "sequence: $sequence"
echo "epsilon: $epsilon"

# Create the target dir for sequence
target_dir="${current_simulation_rg_dir}/b2-1000/${seqid}"
mkdir -p "$target_dir"
cd "$target_dir"

# LAMMPS input file
calcpath="${target_dir}/rg-${seqid}.in"
python3 "${current_simulation_rg_dir}/build_one_chain.py" 2 "$sequence" "polymer_init-${seqid}.data"

# Copy the Rg sim template into target dir, replace epsilon and random seed
cp "${current_simulation_rg_dir}/rg.in" "$calcpath"

# Insert the epsilon value into the LAMMPS input file
sed -i "s/{epsilon}/$epsilon/g" "$calcpath"

random_number=$(( RANDOM % (10000 - 100 + 1) + 100 ))
sed -i "s/^variable seed1 equal .*/variable seed1 equal $random_number/" "$calcpath"
sed -i "s/polymer_init.data/polymer_init-${seqid}.data/" "$calcpath"
sed -i "s/output$/output ${seqid}/" "$calcpath"
sed -i "s/file rg\.dat/file rg-${seqid}\.dat/" "$calcpath"

# Run simulation
echo "Running Rg simulation for $seqid..."
srun $HOME/software/lammps_umbrella/lammps-stable_29Sep2021/build_us/lmp_della -nc -i "$calcpath" -log none > "${seqid}.log" 2>&1
wait
