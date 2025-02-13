#!/bin/bash
#SBATCH --job-name=sample_configs  # Job name
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --cpus-per-task=1         # Number of CPU cores per task
#SBATCH --mem-per-cpu=1500M       # Memory allocation per CPU core
#SBATCH --time=23:59:00           # Time limit (HH:MM:SS)
#SBATCH --array=1-72              # Array of 72 jobs (each for a different sequence)
#SBATCH --constraint=cascade,skylake  # Architecture constraints

# This script samples two-chain polymer configurations from converged ABF simulations using LAMMPS.
# The sampled configurations are used for building contact maps. 
# The key step here is identifying **decorrelated configurations within a specified center-of-mass (COM) distance range**,
# which ensures that the sampled configurations are both relevant (in-target) and statistically independent.
# Two helper scripts are used for this:
# 1. **extract_target_configs.py**: identifies the largest set of decorrelated configurations within a target COM distance range.
# 2. **copy_target_configs.sh**: copies the in-target configurations to a specified output directory.

ulimit -s unlimited
ulimit -c unlimited

module purge
module load intel/2021.1
module load intel-mpi/gcc/2021.1.1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Generalized paths for data and simulation scripts
data_dir="<path-to-data-files>"                     # Directory containing epsilon values and sequence information
epsilon_file="$data_dir/combined_epsilon_400.txt"   # File with seqID and corresponding epsilon values
seqs_file="$data_dir/b2_400_combined_seqs-ps.lst"  # File with seqID and sequence data
simulation_b2_dir="<path-to-simulation-scripts>/b2" # Directory for B2-matching simulation scripts
sample_configs_dir="<path-to-sampling-configs>"     # Base directory containing template input file for LAMMPS
output_dir="$sample_configs_dir/b2-400"             # B2-matched dataset dir to store sampled configs in b2-400/{seqid} subdir

# Create temp scratch dir to save space
ScratchDir="/tmp/${SLURM_JOB_ID}"
mkdir -p ${ScratchDir}

line_number=$((SLURM_ARRAY_TASK_ID))
subdir_line=$(sed -n "${line_number}p" "$seqs_file")

# Extract seqID and sequence from the current line
seqid=$(echo "$subdir_line" | awk '{print $1}')
sequence=$(echo "$subdir_line" | awk '{print $2}')
epsilon=$(awk -v id="$seqid" '$1 == id { print $2 }' "$epsilon_file")

echo "seqID: $seqid"
echo "sequence: $sequence"
echo "epsilon: $epsilon"

seqid_output_dir="${output_dir}/${seqid}"
mkdir -p "${seqid_output_dir}/target-configs/${SLURM_JOB_ID}"

# Copy the LAMMPS input file template, replace epsilon value and random seed
calcpath="${seqid_output_dir}/sample_configs_for_contact_map.in"
cp "${sample_configs_dir}/sample_configs_for_contact_map.in" "$calcpath"

random_number=$(( RANDOM % (10000 - 100 + 1) + 100 ))
sed -i "s/{random_number}/$random_number/" "$calcpath"
sed -i "s/{epsilon}/$epsilon/g" "$calcpath"
sed -i "s/output$/output ${seqid}/" "$calcpath"

if [ ! -f "${seqid_output_dir}/colvars.inp" ]; then
    cp "${simulation_b2_dir}/colvars.inp" "${seqid_output_dir}"
fi

# Copy files to the temp Scratch dir
mkdir -p ${ScratchDir}/$seqid/configs
cp -r "${seqid_output_dir}" ${ScratchDir}
cp "${sample_configs_dir}/copy_target_configs.sh" ${ScratchDir}/$seqid
cp "${sample_configs_dir}/extract_target_configs.py" ${ScratchDir}/$seqid

cd "${ScratchDir}/$seqid"

# Run LAMMPS simulation
echo "Running LAMMPS simulation to sample configs for for $seqid..."
srun $HOME/software/lammps_umbrella/lammps-stable_29Sep2021/build_us/lmp_della -nc -in "$calcpath" -log none > "${seqid}.log" 2>&1
wait

# Run py script to identify the largest set of decorrelated in-target configurations
# This script scans the .colvars.traj file and selects configurations with a COM distance in the specified range (e.g., 1.0 sigma <= COM <= 5.0 sigma)
echo "Identifying decorrelated configurations within the target COM distance range..."
python3 extract_target_configs.py "${seqid}"

# Run the bash script to copy the identified configurations to the output dir
# This script reads the winning_sequence.colvars.traj file (created by extract_target_configs.py) and copies the corresponding .atom files
echo "Copying decorrelated configurations..."
bash copy_target_configs.sh

# List the contents of target-configs/ for verification
ls -lh target-configs/

# Save all in-target decorrelated configurations to a final output dir (named by SLURM_JOB_ID)
#  This avoids file overwriting when running parallel jobs, because dumped files may have the same names across different runs. Later, all configurations will be collected from these subdirectories to build the contact map.

echo "Copying decorrelated configurations to final output directory (${seqid_output_dir}/target-configs/${SLURM_JOB_ID})..."
cp target-configs/*.atom "${seqid_output_dir}/target-configs/${SLURM_JOB_ID}"

echo "Configurations sampled for seqID: $seqid"
