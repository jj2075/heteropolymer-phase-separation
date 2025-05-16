#!/bin/bash

# Source data directories
source_dir="/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation"

# Paths to sequence lists (csv)
seq_list_400="$source_dir/sequences/hp_b2_400_dataset.csv"
seq_list_1000="$source_dir/sequences/hp_b2_1000_dataset.csv"
seq_list_idp="$source_dir/sequences/idp_challenge_dataset.csv"

# Define contact cutoff distances for each model
HP_CUTOFFS=("2.0" "2.5" "3.0" "3.5" "4.0")
IDP_CUTOFFS=("12" "16" "20" "24" "28" "32")

# STEP 1: Consolidate feature numpy arrays for each B2 value from raw data
b2_value_400=400
b2_value_1000=1000

# Output directories for feature arrays
output_dir_400="feature_npy_arrays/hp_b2_400"
output_dir_1000="feature_npy_arrays/hp_b2_1000"
output_dir_idp="feature_npy_arrays/idp_challenge_set"

mkdir -p "$output_dir_400"
mkdir -p "$output_dir_1000"
mkdir -p "$output_dir_idp"

# Call the py script for HP (b2=400)
for cutoff in "${HP_CUTOFFS[@]}"; do
    python3 consolidate_feature_arrays.py \
        --seq_list_path "$seq_list_400" \
        --b2_value "$b2_value_400" \
        --rg_dir "$source_dir/single_chain_rg/hp_b2_400" \
        --mean_dir "$source_dir/contact_map_statistics/mean/hp_b2_400/cutoff_${cutoff}_sigma" \
        --variance_dir "$source_dir/contact_map_statistics/variance/hp_b2_400/cutoff_${cutoff}_sigma" \
        --partial_sums_dir "$source_dir/contact_map_fourier_outputs/partial_power_sums/hp_b2_400/cutoff_${cutoff}_sigma" \
        --output_dir "$output_dir_400/cutoff_${cutoff}" \
        --nyquist_frequency 10 \
        --cutoff_distance "${cutoff}"
done

# Call the py script for HP (b2=1000)
for cutoff in "${HP_CUTOFFS[@]}"; do
    python3 consolidate_feature_arrays.py \
        --seq_list_path "$seq_list_1000" \
        --b2_value "$b2_value_1000" \
        --rg_dir "$source_dir/single_chain_rg/hp_b2_1000" \
        --mean_dir "$source_dir/contact_map_statistics/mean/hp_b2_1000/cutoff_${cutoff}_sigma" \
        --variance_dir "$source_dir/contact_map_statistics/variance/hp_b2_1000/cutoff_${cutoff}_sigma" \
        --partial_sums_dir "$source_dir/contact_map_fourier_outputs/partial_power_sums/hp_b2_1000/cutoff_${cutoff}_sigma" \
        --output_dir "$output_dir_1000/cutoff_${cutoff}" \
        --nyquist_frequency 10 \
        --cutoff_distance "${cutoff}"
done

# Call the py script for the IDP model (per-sequence b2 values from CSV)
for cutoff in "${IDP_CUTOFFS[@]}"; do
    python3 consolidate_feature_arrays.py \
        --seq_list_path "$seq_list_idp" \
        --rg_dir "$source_dir/single_chain_rg/idp_challenge_set" \
        --mean_dir "$source_dir/contact_map_statistics/mean/idp_challenge_set/cutoff_${cutoff}_angstroms" \
        --variance_dir "$source_dir/contact_map_statistics/variance/idp_challenge_set/cutoff_${cutoff}_angstroms" \
        --partial_sums_dir "$source_dir/contact_map_fourier_outputs/partial_power_sums/idp_challenge_set/cutoff_${cutoff}_angstroms" \
        --output_dir "$output_dir_idp/cutoff_${cutoff}" \
        --nyquist_frequency 10 \
        --cutoff_distance "${cutoff}"
done

# STEP 2: Combine feature arrays for B2-matched datasets in order to train model on mixed B2 data
combined_output_dir="feature_npy_arrays/hp_b2_mixed"
mkdir -p "$combined_output_dir"

for cutoff in "${HP_CUTOFFS[@]}"; do
    python3 combine_b2_matched_feature_arrays.py \
        --b2_dir_400 "$output_dir_400/cutoff_${cutoff}" \
        --b2_dir_1000 "$output_dir_1000/cutoff_${cutoff}" \
        --rg_dir_400 "$output_dir_400/cutoff_${cutoff}" \
        --rg_dir_1000 "$output_dir_1000/cutoff_${cutoff}" \
        --mean_dir_400 "$output_dir_400/cutoff_${cutoff}" \
        --mean_dir_1000 "$output_dir_1000/cutoff_${cutoff}" \
        --variance_dir_400 "$output_dir_400/cutoff_${cutoff}" \
        --variance_dir_1000 "$output_dir_1000/cutoff_${cutoff}" \
        --partial_sums_dir_400 "$output_dir_400/cutoff_${cutoff}/partial_powers" \
        --partial_sums_dir_1000 "$output_dir_1000/cutoff_${cutoff}/partial_powers" \
        --output_dir "$combined_output_dir/cutoff_${cutoff}" \
        --nyquist_frequency 10
done

