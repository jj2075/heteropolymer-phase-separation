import os
import csv
import json
import random
import numpy as np
import argparse
from glob import glob
from collections import OrderedDict
import warnings
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.metrics import (
    make_scorer, log_loss, accuracy_score, confusion_matrix, precision_score,
    recall_score, f1_score, roc_auc_score, roc_curve, auc)
warnings.filterwarnings("ignore")

"""
Script for training and evaluating a NON-split-sum logistic regression model. The model can use any combination of the following features: B_22, contact-map mean, contact-map variance, single-chain radius of gyration.
Model training uses a 80/20 train-test split.
"""

def load_combination_mapping(json_file_path):
    """
    load JSON file containing model-feature combination
    """
    with open(json_file_path, 'r') as json_file:
        combination_variable_mapping = json.load(json_file)
    combination_list = list(combination_variable_mapping.items())

    combination_str, combo_name = combination_list[0]
    combination_str = combination_str.strip('()').replace("'", "").replace(" ", "")
    if combination_str.endswith(','):
        combination_str = combination_str[:-1]

    return combination_str, combo_name

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
    # Create a dictionary to map filenames to their full paths
    file_path_map = {}
    for root, _, files in os.walk(feature_data_dir):
        for file in files:
            if file.endswith(".npy"):
                file_path_map[file] = os.path.join(root, file)

    # Match npy files with their full paths
    npy_files = []
    for npy_file in combination_str.split(','):
        npy_file = npy_file.strip()
        if npy_file in file_path_map:
            npy_files.append(file_path_map[npy_file])
        else:
            raise FileNotFoundError(f"The feature file {npy_file} is not found in {feature_data_dir}")

    return npy_files

# get the next available trial index for this combination to avoid overwriting results
def get_next_trial_index(model_results_dir):
    existing_files = glob(os.path.join(model_results_dir, 'roc_data_trial_*.json'))
    if not existing_files:
        return 0  # If no files exist, start with index 0
    existing_indices = [int(os.path.splitext(os.path.basename(f))[0].split('_')[-1]) for f in existing_files]
    return max(existing_indices) + 1  # Return the next index

def load_data(npy_files, seq_ps_file):
    """
    Load feature data from .npy files, match them to sequence labels (1 for phase-separating, 0 for non-phase-separating).
    Parameters:
       - npy_files (list): list of .npy file paths containing structured 2D seqID-feature data
        - seq_ps_file (str): path to the file containing sequence labels

    Returns:
        - ordered_feature_dicts (list): list of ordered dictionaries mapping common seqIDs to feature values
        labels_dict (dict): dict mapping for common sequence IDs to truth labels
    """
    loaded_features = [np.load(npy_file, allow_pickle=True) for npy_file in npy_files]
    feature_dicts = [{seq: feature for seq, feature in loaded_feature} for loaded_feature in loaded_features]

    common_seqs = set.intersection(*[set(feature_dict.keys()) for feature_dict in feature_dicts])
    if not common_seqs:
        raise ValueError("No common seqIDs found between feature arrays and labels")
    ordered_feature_dicts = [OrderedDict(sorted((k, feature_dict[k]) for k in common_seqs)) for feature_dict in feature_dicts]

    labels_dict = {}
    with open(seq_ps_file, 'r') as f:
        reader = csv.reader(f) # for CSV file
        next(reader)
        for row in reader:
            seq_id, label = row[0], int(row[1])  # SeqID is in column 1, Phase-Sep is in column 2
            if seq_id in common_seqs:
                labels_dict[seq_id] = label

    return ordered_feature_dicts, labels_dict

def preprocess_data(feature_dicts, labels_dict, test_size=0.2, random_state=None):
    X = []
    y = []
    indices = []

    for i, key in enumerate(labels_dict.keys()):
        X.append([float(feature_dict[key]) for feature_dict in feature_dicts if key in feature_dict])
        y.append(labels_dict[key])
        indices.append(i)

    X_train, X_test, y_train, y_test, train_indices, test_indices = train_test_split(
        X, y, indices, test_size=test_size, random_state=random_state
    )

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    y_train = np.array(y_train).astype(int)
    y_test = np.array(y_test).astype(int)

    print("Training features shape:", X_train_scaled.shape)
    print("Test features shape:", X_test_scaled.shape)

    return X_train_scaled, X_test_scaled, y_train, y_test, train_indices, test_indices

def train_and_evaluate_model(X_train, X_test, y_train, y_test, train_indices, test_indices):
    """
    Train and evaluate an unregularized logistic regression model.
    """
    model = LogisticRegression(penalty='none', solver='lbfgs')
    model.fit(X_train, y_train)
    y_pred_prob_train = model.predict_proba(X_train)[:, 1]
    y_pred_prob_test = model.predict_proba(X_test)[:, 1]

    y_train_pred = model.predict(X_train).tolist()
    y_test_pred = model.predict(X_test).tolist()

    train_auc = roc_auc_score(y_train, y_pred_prob_train)
    test_auc = roc_auc_score(y_test, y_pred_prob_test)

    evaluation_details = {
        "train_auc": train_auc,
        "test_auc": test_auc,
        "train_indices": train_indices,
        "test_indices": test_indices,
        "coefficients": model.coef_.tolist(),
        "intercept": model.intercept_.tolist(),
        "train_predicted_labels": y_train_pred,
        "true_labels_test": y_test.tolist(),
        "test_predicted_labels": y_test_pred,
        "train_predicted_probabilities": y_pred_prob_train.tolist(),
        "test_predicted_probabilities": y_pred_prob_test.tolist()
    }

    return evaluation_details

def save_model_performance(model_results_dir, evaluation_details, trial_index):
    print("Writing model evaluation details...")
    with open(os.path.join(model_results_dir, f"roc_data_trial_{trial_index}.json"), 'w') as file:
        json.dump(evaluation_details, file)

def run_trial(npy_files, seq_ps_file, model_results_dir, trial_index, random_state):
    os.makedirs(model_results_dir, exist_ok=True)
    print(f"Running trial {trial_index} with random_state {random_state}")

    # Load data
    feature_dicts, labels_dict = load_data(npy_files, seq_ps_file)

    # Split into test and train sets
    X_train, X_test, y_train, y_test, train_indices, test_indices = preprocess_data(feature_dicts, labels_dict, random_state=random_state)

    # Train and evaluate model
    evaluation_details = train_and_evaluate_model(X_train, X_test, y_train, y_test, train_indices, test_indices)
    save_model_performance(model_results_dir, evaluation_details, trial_index)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train and evaluate a non-split-sum logistic regression model")
    parser.add_argument("--model_feature_mapping_file", type=str, required=True, help="Path to JSON file specifying model-feature mapping")
    parser.add_argument("--feature_data_dir", type=str, required=True, help="Directory containing numpy feature data arrays")
    parser.add_argument("--model_results_dir", type=str, required=True, help="Folder to store model performance results")
    parser.add_argument("--seq_ps_file", type=str, required=True, help="Path to the file with sequence IDs and phase-separation labels")
    parser.add_argument("--random_state", type=int, required=True, help="Random state for test-train split")
    args = parser.parse_args()

    # Use provided JSON file path and load its specified model-feature mapping
    model_feature_mapping_file = args.model_feature_mapping_file
    combination_str, combo_name = load_combination_mapping(model_feature_mapping_file)
    npy_files = get_npy_files(combination_str, args.feature_data_dir)

    # Get next trial index (for writing results), then train and evaluate the model
    trial_index = get_next_trial_index(args.model_results_dir)
    run_trial(npy_files, args.seq_ps_file, args.model_results_dir, trial_index, args.random_state)
