#!/bin/bash
#SBATCH --job-name=logr-model      # job name
#SBATCH --nodes=1                  # number of nodes
#SBATCH --ntasks=1                 # number of tasks
#SBATCH --cpus-per-task=1          # cpu cores per task
#SBATCH --mem-per-cpu=200M         # memory per cpu core
#SBATCH --time=00:30:00            # job run time limit (hh:mm:ss)

# IDP contact distance cutoff used in the paper (24 Angstroms)
cutoff_distance=24

# Paths for inputs and outputs
model_feature_mapping_file="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_feature_mappings/b2_rg_mean_variance.json"
feature_data_dir="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/data_processing/feature_npy_arrays/idp_challenge_set/cutoff_${cutoff_distance}/"
seq_ps_file="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/idp_challenge_set_phase_separation_labels.csv"
model_results_dir="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_results/b2-rg-mean-variance/idp_challenge_set/cutoff_${cutoff_distance}/"
mkdir -p "$model_results_dir"

# Loop for 1 trials
for i in {1..184}; do
    # Generate a random number between 0 and 10000 for random state for test-train split
    random_state=$((RANDOM % 10001))

    # Train and evaluate the model
    python3 train_test_any_non_split_sum_model.py \
        --model_feature_mapping_file "$model_feature_mapping_file" \
        --feature_data_dir "$feature_data_dir" \
        --seq_ps_file "$seq_ps_file" \
        --model_results_dir "$model_results_dir" \
        --random_state "$random_state"
done
