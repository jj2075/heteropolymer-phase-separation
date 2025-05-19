import os
import json

# Define base directories/files for HP (b2-mixed) model and IDP model
base_dir = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_results"
seq_ps_file_idp = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/idp_challenge_set_phase_separation_labels.csv"

commands = []

# IDP b2-rg-mean-var model, split-sum model has same performance)
idp_cutoffs = ["16", "20", "24", "28", "32"]
for cutoff in idp_cutoffs:
    cmd = f"""python3 analyze_logistic_regression_model_results.py \\
  --model_results_dir "{base_dir}/b2-rg-mean-variance/idp_challenge_set/cutoff_{cutoff}" \\
  --seq_ps_file "{seq_ps_file_idp}" \\
  --save_auc"""
    commands.append(cmd)

import subprocess
for cmd in commands:
    print(f"Running:\n{cmd}\n{'-'*60}")
    subprocess.run(cmd, shell=True, check=True)
