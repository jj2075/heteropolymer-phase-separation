#!/bin/bash
#SBATCH --job-name=split-sum-model  # job name
#SBATCH --nodes=1                   # number of nodes
#SBATCH --ntasks=1                  # number of tasks
#SBATCH --cpus-per-task=1           # cpu cores per task
#SBATCH --mem-per-cpu=200M          # memory per cpu core
#SBATCH --time=01:05:00             # job run time limit (hh:mm:ss)
#SBATCH --array=0-64                # 65 split-sum models to run

# Contact distance cutoff used in the paper (3.0 sigma)
cutoff_distance=3.0
cutoff_distance=2.5
cutoff_distance=3.5
cutoff_distance=4.0

# Gamma value for this model run
gamma=6.5e-05

# Paths for inputs and outputs
model_feature_mapping_file="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_feature_mappings/split_sums_b2_mixed.json"
feature_data_dir="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/data_processing/feature_npy_arrays/hp_b2_mixed/cutoff_${cutoff_distance}/"
seq_ps_file="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/hp_b2_mixed_phase_separation_labels.csv"
model_results_dir="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_results/split_sum_models/hp_b2_mixed/cutoff_${cutoff_distance}/"
mkdir -p "$model_results_dir"

# Loop for 1 trials
for i in {1..185}; do
    # Generate a random number between 0 and 10000 for random_state
    random_state=$((RANDOM % 10001))

    # Run the model
    python3 train_test_a_split_sum_model.py \
        --entry_index "$((SLURM_ARRAY_TASK_ID))" \
        --model_feature_mapping_file "$model_feature_mapping_file" \
        --feature_data_dir "$feature_data_dir" \
        --seq_ps_file "$seq_ps_file" \
        --model_results_dir "$model_results_dir" \
        --random_state "$random_state" \
        --gamma "$gamma"
done
