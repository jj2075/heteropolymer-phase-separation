import os
import json

# Define base directories/files for HP (b2-mixed) model and IDP model
base_dir = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_results"
mapping_file_hp = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_feature_mappings/split_sums_b2_mixed.json"
seq_ps_file_hp = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/hp_b2_mixed_phase_separation_labels.csv"
seq_ps_file_idp = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/idp_challenge_set_phase_separation_labels.csv"

# Load split-sum model names e.g., power-split-0-4
with open(mapping_file_hp, 'r') as f:
    combo_dict_hp = json.load(f)
combo_names_hp = list(combo_dict_hp.values())

commands = []

# HP B2-mixed model - split-sum model (best)
hp_cutoffs = ["2.5", "3.0", "3.5", "4.0"]
for cutoff in hp_cutoffs:
    for combo in combo_names_hp:
        cmd = f"""python3 analyze_logistic_regression_model_results.py \\
  --model_results_dir "{base_dir}/split_sum_models/hp_b2_mixed/cutoff_{cutoff}" \\
  --seq_ps_file "{seq_ps_file_hp}" \\
  --split_sum \\
  --combo_name "{combo}" \\
  --save_auc"""
        commands.append(cmd)

# IDP b2-rg-mean-var model
idp_cutoffs = ["16", "20", "24", "28", "32"]
for cutoff in idp_cutoffs:
    cmd = f"""python3 analyze_logistic_regression_model_results.py \\
  --model_results_dir "{base_dir}/b2-rg-mean-variance/idp_challenge_set/cutoff_{cutoff}" \\
  --seq_ps_file "{seq_ps_file_idp}" \\
  --save_auc"""
    commands.append(cmd)

# IDP split-sum model
idp_cutoffs = ["20", "24", "28"]
for cutoff in idp_cutoffs:
    for combo in combo_names_hp:
        cmd = f"""python3 analyze_logistic_regression_model_results.py \\
  --model_results_dir "{base_dir}/split_sum_models/idp_challenge_set/cutoff_{cutoff}" \\
  --seq_ps_file "{seq_ps_file_idp}" \\
  --split_sum \\
  --combo_name "{combo}" \\
  --save_auc"""
        commands.append(cmd)


import subprocess
for cmd in commands:
    print(f"Running:\n{cmd}\n{'-'*60}")
    subprocess.run(cmd, shell=True, check=True)
