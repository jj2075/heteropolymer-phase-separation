import os
import json

# Define base paths
base_dir = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_results"
hp_cutoffs = ["2.5", "3.0", "3.5", "4.0"]
idp_cutoffs = ["16", "20", "24", "28", "32"]
model_feature_mapping_file = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation/logistic_regression/model_feature_mappings/split_sums_b2_mixed.json"

hp_summary = []
idp_summary = []

# Load JSON mapping
with open(model_feature_mapping_file, "r") as f:
    model_mappings = json.load(f)

dir_names = set(model_mappings.values())

# Process HP results
for cutoff in hp_cutoffs:
    cutoff_dir = os.path.join(base_dir, "split_sum_models", "hp_b2_mixed", f"cutoff_{cutoff}")
    best_auc = -1
    best_stderr = float("inf")
    best_model = None

    if not os.path.exists(cutoff_dir):
        continue

    for model in os.listdir(cutoff_dir):
        if model not in dir_names:
            continue

        path = os.path.join(cutoff_dir, model, "avg_model_test_auc.dat")
        if os.path.exists(path):
            with open(path, "r") as f:
                line = f.readline().strip()
                if "±" in line or "+=" in line:
                    try:
                        if "±" in line:
                            auc_str, stderr_str = line.split("±")
                        elif "+=" in line:
                            auc_str, stderr_str = line.split("+=")
                        auc = float(auc_str.strip())
                        stderr = float(stderr_str.strip())
                        if auc > best_auc or (auc == best_auc and stderr < best_stderr):
                            best_auc = auc
                            best_stderr = stderr
                            best_model = model
                    except ValueError:
                        print(f"Error parsing file {path}: {line}")
                        continue
    if best_model:
        hp_summary.append((cutoff, best_model, best_auc, best_stderr))

# Process IDP results
for cutoff in idp_cutoffs:
    model_dir = os.path.join(base_dir, "b2-rg-mean-variance", "idp_challenge_set", f"cutoff_{cutoff}")
    path = os.path.join(model_dir, "avg_model_test_auc.dat")
    best_auc = -1
    best_stderr = float("inf")

    if os.path.exists(path):
        with open(path, "r") as f:
            line = f.readline().strip()
            # Handle both `±` and `+=` formats
            if "±" in line or "+=" in line:
                try:
                    if "±" in line:
                        auc_str, stderr_str = line.split("±")
                    elif "+=" in line:
                        auc_str, stderr_str = line.split("+=")
                    auc = float(auc_str.strip())
                    stderr = float(stderr_str.strip())
                    if auc > best_auc or (auc == best_auc and stderr < best_stderr):
                        best_auc = auc
                        best_stderr = stderr
                except ValueError:
                    print(f"Error parsing file {path}: {line}")
                    continue
    if best_auc > -1:
        idp_summary.append((cutoff, best_auc, best_stderr))

if not hp_summary:
    raise ValueError("No HP results were collected. Verify your directories and JSON mappings.")
if not idp_summary:
    raise ValueError("No IDP results were collected. Verify your directories and cutoff values.")

# Print HP summary
print("\nHP Summary Results:")
print(f"{'Cutoff':<10} {'Best Model':<30} {'AUC':<10} {'Std Err':<10}")
for cutoff, model, auc, stderr in hp_summary:
    print(f"{cutoff:<10} {model:<30} {auc:<10.4f} {stderr:<10.4f}")

# Print IDP summary
print("\nIDP Summary Results:")
print(f"{'Cutoff':<10} {'AUC':<10} {'Std Err':<10}")
for cutoff, auc, stderr in idp_summary:
    print(f"{cutoff:<10} {auc:<10.4f} {stderr:<10.4f}")
