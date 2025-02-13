import numpy as np
import os
import argparse
import random
import json
from glob import glob
from collections import OrderedDict
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.metrics import (
    make_scorer, log_loss, accuracy_score, confusion_matrix, precision_score,
    recall_score, f1_score, roc_auc_score, roc_curve, auc)
)
import warnings
warnings.filterwarnings("ignore")

"""
Script for training and evaluating split-sum logistic regression models for the B2-combined dataset.
Models are trained using an 80/20 train-test split with features representing partial sums of the power spectrum (A, B)
and the second virial coefficient, $B_{22}$(F). Baseline performance is established using the simpler (C, F) model,
where C represents the total variance (sum of A and B).

Model selection is based on maximizing the training AUC, with a penalty applied to discourage models
with marginal training improvements that may not generalize well, i.e., avoid overfitting.
The best-performing model on the training set is evaluated on the test set.
"""

def load_combination_mapping(json_file_path, entry_index):
    """
    Reads the JSON file, extracts the (entry_index)-th combination.
    The JSON file "combination_split_sums_powers.json" maps feature combinations to combination names.
    Each key is a tuple of three .npy file paths:
     - partial sum of the power spectrum up to a certain wavenumber (e.g., "sum_up_to_include_0_1_power.npy")
     - partial sum of the power spectrum for modes above that wavenumber (e.g., "sum_above_0_1_power.npy")
     - the B2 value (e.g., "b2.npy")
    Value is the combination name (str) (e.g., "power-split-0-1").

    Returns: combination_str (str) representing the tuple of .npy file paths and combo_name (str)
    representing name associated with the feature combination (e.g., "power-split-0-1").
    """
    with open(json_file_path, 'r') as json_file:
        combination_variable_mapping = json.load(json_file)
    combination_list = list(combination_variable_mapping.items())

    combination_str, combo_name = combination_list[entry_index]
    combination_str = combination_str.strip('()').replace("'", "").replace(" ", "")
    if combination_str.endswith(','):
        combination_str = combination_str[:-1]

    return combination_str, combo_name

def get_npy_files(combination_str, feature_data_dir):

def get_npy_files(combination_str, feature_data_dir):
    """
    Given a specified combination of features, load the corresponding .npy file paths.
    Parameters:
        - combination_str (str): comma-separated list of .npy filenames (e.g., "feature1.npy, feature2.npy, feature3.npy").
        - feature_data_dir (str): base directory containing feature data
    Returns:
       - npy_files (list): List of file paths for the .npy files in feature_data_dir

    Each .npy file is a structured 2D array: row holds seqIDs, col holds feature values.
    """
    npy_files = []
    for npy_file in combination_str.split(','):
        npy_file = npy_file.strip()
        npy_file_path = os.path.join(feature_data_dir, npy_file)
        
        if os.path.exists(npy_file_path):
            npy_files.append(npy_file_path)
        else:
            print(f"Warning: {npy_file_path} not found. Skipping.")
    
    return npy_files

# get the next available trial index for this combination to avoid overwriting results
def get_next_trial_index(results_dir):
    existing_files = glob(os.path.join(results_dir, 'roc_data_trial_*.json'))
    if not existing_files:
        return 0  # If no files exist, start with index 0
    existing_indices = [int(os.path.splitext(os.path.basename(f))[0].split('_')[-1]) for f in existing_files]
    return max(existing_indices) + 1  # Return the next index

def load_data(npy_files, seq_ps_file):
    """
    Load feature data from .npy files, match them to sequence labels (1 for phase-separating, 0 for aggregating).
    Parameters:
       - npy_files (list): list of .npy file paths containing feature data
        - seq_ps_file (str): path to the file containing sequence labels

    Returns:
        - ordered_feature_dicts (list): list of ordered dictionaries mapping common seqIDs to feature values
        labels_dict (dict): dict mapping for common sequence IDs to truth labels
    """
    loaded_features = [np.load(npy_file, allow_pickle=True) for npy_file in npy_files]
    feature_dicts = [{seq: feature for seq, feature in loaded_feature} for loaded_feature in loaded_features]

    common_seqs = set.intersection(*[set(feature_dict.keys()) for feature_dict in feature_dicts])
    ordered_feature_dicts = [OrderedDict(sorted((k, feature_dict[k]) for k in common_seqs)) for feature_dict in feature_dicts]

    labels_dict = {}
    with open(seq_ps_file, 'r') as f:
        for line in f:
            cols = line.split()
            seq, label = cols[0], int(cols[2])
            if seq in common_seqs:
                labels_dict[seq] = label - 1

    return ordered_feature_dicts, labels_dict

def get_roc_details(model_name, model, X_train, y_train, X_test, y_test):
    """
    Generate ROC details for a given model.
    Params: model_name, the log-reg model, X_train, y_train, X_test, y_test
    Returns: dictionary of AUC scores, model coefficients, intercept, predicted labels, and probabilities. 
    """
    y_train_pred_prob = model.predict_proba(X_train)[:, 1]
    y_test_pred_prob = model.predict_proba(X_test)[:, 1]

    y_train_pred = model.predict(X_train).tolist()
    y_test_pred = model.predict(X_test).tolist()

    train_auc = roc_auc_score(y_train, y_train_pred_prob)
    test_auc = roc_auc_score(y_test, y_test_pred_prob)

    roc_details = {
        "model_name": model_name,
        "train_auc": train_auc,
        "test_auc": test_auc,
        "coefficients": model.coef_.tolist(),
        "intercept": model.intercept_.tolist(),
        "train_predicted_labels": y_train_pred,
        "test_predicted_labels": y_test_pred,
        "train_predicted_probabilities": y_train_pred_prob.tolist(),
        "test_predicted_probabilities": y_test_pred_prob.tolist()
    }

    return roc_details


def build_cf_data(X_train_raw, X_test_raw):
    """
    Build "C-F" data by summing features A and B before applying standardization (StandardScaler).
    Parameters:
      - X_train_raw (array): training features (A, B, F), shape (n_samples, 3)
      - X_test_raw (array): test features (A, B, F), shape (n_samples, 3)
    Returns:
      - X_train_cf_raw (array): training features for C-F model (C = A + B, F)
      - X_test_cf_raw (array): test features for C-F model (C = A + B, F)
    """
    AplusB_train = (X_train_raw[:, 0] + X_train_raw[:, 1]).reshape(-1, 1)
    AplusB_test = (X_test_raw[:, 0] + X_test_raw[:, 1]).reshape(-1, 1)
    F_train_raw = X_train_raw[:, 2].reshape(-1, 1)
    F_test_raw = X_test_raw[:, 2].reshape(-1, 1)

    X_train_cf_raw = np.hstack([AplusB_train, F_train_raw])
    X_test_cf_raw = np.hstack([AplusB_test, F_test_raw])

    return X_train_cf_raw, X_test_cf_raw


def rescale_cf_to_abf(scaler_cf, scaler_abf, coef_cf, intercept_cf):
    """
    Rescale C-F model coefficients to derive equivalent AB-F model coefficients.
    Parameters:
      - scaler_cf (StandardScaler): scaler used for C-F standardization
      - scaler_abf (StandardScaler): scaler used for AB-F standardization
      - coef_cf (array): coefficients of the C-F model
      - intercept_cf (float): intercept of the C-F model
    Returns:
      - derived_coefs (array): coefficients for the AB-F model (alpha_a, alpha_b, alpha_f)
      - intercept_abf_derived (float): intercept for the AB-F model
    """
    # C-F means and scales
    mean_c = scaler_cf.mean_[0]  # mean(A+B)
    mean_fcf = scaler_cf.mean_[1]  # mean(F) used in C-F scaling
    std_c = scaler_cf.scale_[0]  # std(A+B)
    std_fcf = scaler_cf.scale_[1]  # std(F)
    beta_c = coef_cf[0]  # coefficient for scaled(A+B)
    beta_f = coef_cf[1]  # coefficient for scaled(F)

    # AB-F means and scales
    mean_a = scaler_abf.mean_[0]
    mean_b = scaler_abf.mean_[1]
    mean_fabf = scaler_abf.mean_[2]
    std_a = scaler_abf.scale_[0]
    std_b = scaler_abf.scale_[1]
    std_fabf = scaler_abf.scale_[2]

    # Compute derived coefficients
    alpha_a = beta_c * (std_a / std_c)
    alpha_b = beta_c * (std_b / std_c)
    alpha_f = beta_f * (std_fabf / std_fcf)

    # Compute derived intercept
    intercept_abf_derived = (
        intercept_cf
        - beta_c * (mean_c / std_c)  # C-F contribution from raw means
        - beta_f * (mean_fcf / std_fcf)
        + alpha_a * (mean_a / std_a)  # AB-F contribution
        + alpha_b * (mean_b / std_b)
        + alpha_f * (mean_fabf / std_fabf)
    )

    return np.array([alpha_a, alpha_b, alpha_f]), intercept_abf_derived

def train_and_evaluate_models(X_train_raw, X_test_raw, y_train, y_test):
    """
    Train and evaluate log regression models for both AB-F and C-F models.
    Steps:
      1. build C-F data by summing features A and B
      2. scale both AB-F and C-F data using StandardScaler
      3. train unregularized log regression models for AB-F and C-F
      4. derive the equivalent AB-F model from the C-F model by rescaling coefficients
    Returns a dict containing:
       - trained AB-F and C-F models
       - derived AB-F model (from C-F)
       - training and test AUCs for each model
       - scaled feature matrices (AB-F and C-F)
    """
    # Step 1: Build C-F data
    X_train_cf_raw, X_test_cf_raw = build_cf_data(X_train_raw, X_test_raw)

    # Step 2: Scale AB-F data
    scaler_abf = StandardScaler()
    X_train_abf = scaler_abf.fit_transform(X_train_raw)
    X_test_abf = scaler_abf.transform(X_test_raw)

    # Step 3: Scale C-F data
    scaler_cf = StandardScaler()
    X_train_cf = scaler_cf.fit_transform(X_train_cf_raw)
    X_test_cf = scaler_cf.transform(X_test_cf_raw)

    # Step 4: Train AB-F model
    model_abf = LogisticRegression(penalty='none', solver='lbfgs')
    model_abf.fit(X_train_abf, y_train)
    auc_abf_train = roc_auc_score(y_train, model_abf.predict_proba(X_train_abf)[:, 1])
    auc_abf_test = roc_auc_score(y_test, model_abf.predict_proba(X_test_abf)[:, 1])

    # Step 5: Train C-F model
    model_cf = LogisticRegression(penalty='none', solver='lbfgs')
    model_cf.fit(X_train_cf, y_train)
    auc_cf_train = roc_auc_score(y_train, model_cf.predict_proba(X_train_cf)[:, 1])
    auc_cf_test = roc_auc_score(y_test, model_cf.predict_proba(X_test_cf)[:, 1])

    # Step 6: Derive AB-F coefficients from C-F
    coef_cf = model_cf.coef_[0]
    intercept_cf = model_cf.intercept_[0]
    coef_abf_derived, intercept_abf_derived = rescale_cf_to_abf(scaler_cf, scaler_abf, coef_cf, intercept_cf)

    # Derived AB-F model
    derived_model_abf = LogisticRegression(penalty='none', solver='lbfgs')
    derived_model_abf.classes_ = np.array([0, 1])
    derived_model_abf.coef_ = coef_abf_derived.reshape(1, -1)
    derived_model_abf.intercept_ = np.array([intercept_abf_derived])

    # Return models and results
    return {
        "model_abf": model_abf,
        "model_cf": model_cf,
        "derived_model_abf": derived_model_abf,
        "auc_abf_train": auc_abf_train,
        "auc_abf_test": auc_abf_test,
        "auc_cf_train": auc_cf_train,
        "auc_cf_test": auc_cf_test,
        "X_train_abf": X_train_abf,
        "X_test_abf": X_test_abf,
        "X_train_cf": X_train_cf,
        "X_test_cf": X_test_cf
    }

def compute_penalized_auc(auc_train_ab, auc_train_c, gamma=0.000065):
    """
    Compute penalized AUC for model selection.
    Applies a penalty to AUC_train_ab if the difference between AUCs is small, favoring model C.
    """
    epsilon = 1e-5
    lambda_penalty = gamma / (abs(auc_train_ab - auc_train_c) + epsilon)
    penalized_auc = auc_train_ab - lambda_penalty

    return penalized_auc

def run_trial(npy_files, seq_ps_file, specific_results_dir, trial_index, random_state, gamma_factor):
    """
    Train and test split-sum models.
    """
    os.makedirs(specific_results_dir, exist_ok=True)
    print(f"Running trial {trial_index} with random_state {random_state}")

    # Load data
    feature_dicts, labels_dict = load_data(npy_files, seq_ps_file)

    X = []
    y = []
    indices = []

    for i, key in enumerate(labels_dict.keys()):
        # For each sample, extract feature values A, B, and F.
        # These correspond to the partial sums of the power spectrum (A and B) and $B_22$ (F).
        A = float(feature_dicts[0][key]) if key in feature_dicts[0] else 0.0
        B = float(feature_dicts[1][key]) if key in feature_dicts[1] else 0.0
        C = float(feature_dicts[2][key]) if key in feature_dicts[2] else 0.0
        # Append the feature set [A, B, F] to X and the corresponding label to y
        X.append([A, B, C])
        y.append(labels_dict[key])
        indices.append(i)

    # Split X, y into test and train
    X_train_raw, X_test_raw, y_train, y_test, train_indices, test_indices = train_test_split(
        X, y, indices, test_size=0.2, random_state=random_state
    )

    X_train_raw = np.array(X_train_raw)
    X_test_raw = np.array(X_test_raw)
    y_train = np.array(y_train).astype(int)
    y_test = np.array(y_test).astype(int)

    # Train and evaluate models
    results = train_and_evaluate_models(X_train_raw, X_test_raw, y_train, y_test)

    print("[ABF model] AUC train:", results["auc_abf_train"], "AUC test:", results["auc_abf_test"])
    print("[CF model] AUC train:", results["auc_cf_train"], "AUC test:", results["auc_cf_test"])

    # Apply model-complexity penalization to the AB-F model's training AUC
    penalized_auc_abf = compute_penalized_auc(results["auc_abf_train"], results["auc_cf_train"], gamma_factor) # auc_cf_train = derived train

    # Compute ROC curve details for each model: AB-F (split-sum), C-F (simpler baseline), and the derived AB-F model.
    # The derived AB-F model uses scaled coefficients and intercept adjustments to match the predictions of the C-F model.
    roc_abf = get_roc_details("Model_ABF", results["model_abf"], results["X_train_abf"], y_train, results["X_test_abf"], y_test)
    roc_cf = get_roc_details("Model_CF", results["model_cf"], results["X_train_cf"], y_train, results["X_test_cf"], y_test)
    roc_derived = get_roc_details("Derived_ABF", results["derived_model_abf"], results["X_train_abf"], y_train, results["X_test_abf"], y_test)

    all_roc_details = {
        "train_indices": train_indices,
        "test_indices": test_indices,
        "true_labels_test": y_test.tolist(),
        "model_details": [
            roc_abf,
            roc_cf,
            roc_derived
        ]
    }

    roc_abf["penalized_auc_train"] = penalized_auc_abf

    # Pick best model based on penalized training AUC
    best_model = "Model_ABF" if penalized_auc_abf > results["auc_cf_train"] else "Derived_ABF"
    all_roc_details["best_training_model"] = best_model
    all_roc_details["best_test_auc"] = roc_derived["test_auc"] if best_model == "Derived_ABF" else results["auc_abf_test"]

    output_file = f"{specific_results_dir}/roc_data_trial_{trial_index}.json"
    with open(output_file, "w") as f:
        json.dump(all_roc_details, f, indent=2)

    summary = f"""
    random state: {random_state}
    gamma_factor: {gamma_factor}
    --- Model ABF ---
    Coefficients: {results['model_abf'].coef_.ravel().tolist()}
    Intercept: {results['model_abf'].intercept_[0]}
    Train AUC: {results['auc_abf_train']}
    Penalized Train AUC: {penalized_auc_abf}
    Test AUC: {results['auc_abf_test']}

    --- Model CF ---
    Coefficients: {results['model_cf'].coef_.ravel().tolist()}
    Intercept: {results['model_cf'].intercept_[0]}
    Train AUC: {results['auc_cf_train']}
    Test AUC:  {results['auc_cf_test']}

    --- Derived ABF from CF ---
    Coefficients: {results['derived_model_abf'].coef_.ravel().tolist()}
    Intercept: {results['derived_model_abf'].intercept_[0]}
    Train AUC: {roc_derived["train_auc"]}
    Test AUC:  {roc_derived["test_auc"]}

    Best Model: {best_model}
    Best Test AUC: {all_roc_details["best_test_auc"]}
    """
    summary_file = f"{specific_results_dir}/summary_trial_{trial_index}.txt"
    with open(summary_file, "w") as f:
        f.write(summary)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train and evaluate split-sum models")
    parser.add_argument("--entry_index", type=int, required=True, help="idx of JSON file combination, corresponds to (kx,ky) mode")
    parser.add_argument("--features_json_dir", type=str, required=True, help="directory of feature JSON files")
    parser.add_argument("--feature_data_dir", type=str, required=True, help="directory with feature data arrays")
    parser.add_argument("--results_dir", type=str, required=True, help="directory to store test results")
    parser.add_argument("--seq_ps_file", type=str, required=True, help="path to the sequence file")
    parser.add_argument("--random_state", type=int, required=True, help="random state for trial")
    parser.add_argument("--gamma_factor", type=float, default=0.000065, help="Penalization factor for AUC computation")

    args = parser.parse_args()

    json_file_path = os.path.join(args.features_json_dir, "combination_split_sums_powers.json")

    # Load split-sum combination and get npy file paths
    combination_str, combo_name = load_combination_mapping(json_file_path, args.entry_index)
    npy_files = get_npy_files(combination_str, args.feature_data_dir)

    specific_results_dir = os.path.join(args.results_dir, combo_name)

    # get next trial index and run trial
    trial_index = get_next_trial_index(specific_results_dir)
    run_trial(npy_files, args.seq_ps_file, specific_results_dir, trial_index, args.random_state, args.gamma_factor)
