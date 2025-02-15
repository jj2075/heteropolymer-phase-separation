#!/bin/bash
#SBATCH --job-name=b2_400_sim        # Job name
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1                   # Single task
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem-per-cpu=512M           # Memory per CPU core
#SBATCH --time=23:59:00              # Time limit (HH:MM:SS)
#SBATCH --array=1-75                 # Array of 72 jobs (each for a different sequence)
#SBATCH --constraint=cascade,skylake # Architecture constraint

# This script runs adaptive biasing force (ABF) simulations using LAMMPS to compute the potential of mean force (PMF) for B2-matched heteropolymers.
# Each simulation is initialized with random seeds and sequence-specific epsilon values to converge PMFs.

ulimit -s unlimited
ulimit -c unlimited

module purge
module load intel/2021.1
module load intel-mpi/gcc/2021.1.1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Generalized paths
data_dir="<path_to_data_files>"        # Dir containing data files for sequences and epsilon values
simulation_b2_dir="<path_to_simulation_scripts>"  # Dir with LAMMPS simulation scripts and necessary data files
rg_dir="${data_dir}/b2-400/feature-data/rg"       # Dir for pre-computed Rg data

# File with epsilon values (H-H interaction strength) for B2-matched sequences (2nd col: epsilon value)
epsilon_file="${data_dir}/epsilon.csv"

# List of sequences
seqs_file="${data_dir}/b2_400_combined_seqs-ps.lst"

trial_number=$1
line_number=$((SLURM_ARRAY_TASK_ID))

# Get seqID, sequence, and epsilon for the current job from CSV file
subdir_line=$(sed -n "${line_number}p" "$epsilon_file")
seqid=$(echo "$subdir_line" | awk -F ',' '{print $1}')     # seqid in the first column
sequence=$(echo "$subdir_line" | awk -F ',' '{print $2}')  # sequence in the second column
epsilon=$(echo "$subdir_line" | awk -F ',' '{print $4}')   # epsilon value in the fourth column

echo "Simulation for seqID: $seqid"
echo "Sequence: $sequence"
echo "Epsilon: $epsilon"

# Create trial dir
trial_dir="${simulation_b2_dir}/b2-400/${seqid}/trial${trial_number}"
mkdir -p "$trial_dir"
cd "$trial_dir"

# Generate input files for LAMMPS
calcpath="${trial_dir}/calcb2-${seqid}.in"
python3 "${simulation_b2_dir}/build_init.py" 2 "$sequence" "$sequence" "polymer_init-${seqid}.data"

# Copy template files into the trial dir and adjust them
cp "${simulation_b2_dir}/calcb2.in" "$calcpath"
cp "${simulation_b2_dir}/colvars.inp" "${trial_dir}"

# Replace epsilon and random seed in input file
sed -i "s/{epsilon}/$epsilon/g" "$calcpath"
random_number=$(( RANDOM % 9999 + 100 ))
sed -i "s/{random_number}/$random_number/" "$calcpath"
sed -i "s/polymer_init.data/polymer_init-${seqid}.data/" "$calcpath"

pmf_folder="${simulation_b2_dir}/b2-400/pmf"
mkdir -p "$pmf_folder"
pmfpath="${pmf_folder}/${seqid}_trial${trial_number}.pmf"

echo "Running ABF simulation..."
srun <path_to_lammps>/lmp_executable -nc -in "$calcpath" -log none > "${seqid}.log" 2>&1
wait

# Move PMF file to the pmf directory
mv "${seqid}.pmf" "$pmfpath"

# Calculate B2 from the PMF and Rg
# 2*Rg is used to approximate the interaction range between two chains
# in the PMF, defining the scale for a Gaussian envelope that smooths noise in PMF data.
echo "Calculating B2..."
rg1=$(cat "${rg_dir}/avg_rg-${seqid}.dat")
rg2=$rg1
rg_sum=$(echo "$rg1 + $rg2" | bc)
B2=$(python3 "${simulation_b2_dir}/get_pmf.py" "$pmfpath" --Rg_sum "$rg_sum" --seq1id "$seqid" --seq2id "$seqid" | grep B2 | cut -d ' ' -f 4)
echo "$B2" >> "b2-${seqid}-${seqid}.txt"

# Calculate center-of-mass (COM) distance decorrelation time from the colvars trajectory
echo "Calculating decorrelation time..."
traj_file=$(ls *.colvars.traj 2> /dev/null)
if [[ -n "$traj_file" ]]; then
  tail -n 500 "$traj_file" > last_500_traj.txt
fi
python3 "${simulation_b2_dir}/decorrelation.py" "$seqid"
