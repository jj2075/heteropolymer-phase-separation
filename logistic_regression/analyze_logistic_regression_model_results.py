import argparse
import json
import os
import numpy as np
from statistics import mean, stdev
from math import sqrt
import csv

def get_test_split_id(seq_ps_file, sequenceID):
    with open(seq_ps_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        for i, row in enumerate(reader):
            if row[0] == sequenceID:
                return i
    raise ValueError(f"Sequence ID {sequenceID} not found in {seq_ps_file}")

def analyze_non_split_sum(model_dir, test_split_id, save_auc=False):
    aucs = []
    correct_preds = []
    true_labels = []
    predicted_labels = []
    found_in_test = 0

    for i in range(185):
        path = os.path.join(model_dir, f"roc_data_trial_{i}.json")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing trial file: {path}")
        with open(path, 'r') as f:
            data = json.load(f)
            aucs.append(data['test_auc'])

            if test_split_id is not None and test_split_id in data['test_indices']:
                found_in_test += 1
                idx = data['test_indices'].index(test_split_id)
                true = data['true_labels_test'][idx]
                pred = data['test_predicted_labels'][idx]
                correct_preds.append(int(true == pred))
                true_labels.append(true)
                predicted_labels.append(pred)

    avg_auc = round(mean(aucs), 3)
    stderr_auc = round(stdev(aucs) / sqrt(len(aucs)), 3)

    if save_auc:
        out_path = os.path.join(model_dir, "avg_model_test_auc.dat")
        with open(out_path, 'w') as f:
            f.write(f"{avg_auc} += {stderr_auc}\n")
        print(f"saved AUC summary results to: {out_path}")

    if test_split_id is not None:
        print(f"Number of test sets containing this sequence: {found_in_test}")
    print("====================")
    print(f"Average test AUC for phase separation prediction: {avg_auc} ± {stderr_auc}")

    if test_split_id is not None and correct_preds:
        true_label = "Phase separates" if true_labels[0] == 1 else "Does not phase separate"
        avg_accuracy = round(mean(correct_preds), 2)
        stderr_accuracy = round(stdev(correct_preds) / sqrt(len(correct_preds)), 3) if len(correct_preds) > 1 else 0.0
        print(f"\nSequence test/train index = {test_split_id}")
        print(f"Ground truth label: {true_label}")
        print(f"Model accuracy for this sequence: {avg_accuracy} ± {stderr_accuracy}")
        print(f"True labels across test sets: {true_labels}")
        print(f"Predicted labels across test sets: {predicted_labels}")
    elif test_split_id is not None:
        print(f"\nSequence test/train index = {test_split_id} was never in any test set")

    print("====================")

def analyze_split_sum(base_model_dir, combo_name, test_split_id, save_auc=False):
    model_dir = os.path.join(base_model_dir, combo_name)
    aucs = []
    correct_preds = []
    true_labels = []
    found_in_test = 0
    predicted_labels = []

    for i in range(185):
        path = os.path.join(model_dir, f"roc_data_trial_{i}.json")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing trial file: {path}")

        with open(path, 'r') as f:
            data = json.load(f)
            aucs.append(data['best_test_auc'])

            if test_split_id is not None and test_split_id in data['test_indices']:
                found_in_test += 1
                idx = data['test_indices'].index(test_split_id)
                true = data['true_labels_test'][idx]
                true_labels.append(true)

                best_model = data['best_training_model']
                test_preds = None
                for model in data['model_details']:
                    if model['model_name'] == best_model:
                        test_preds = model['test_predicted_labels']
                        break

                if test_preds is None:
                    raise ValueError(f"Best model {best_model} not found in trial {i}")

                predicted = test_preds[idx]
                predicted_labels.append(predicted)
                correct_preds.append(int(predicted == true))

    avg_auc = round(mean(aucs), 2)
    stderr_auc = round(stdev(aucs) / sqrt(len(aucs)), 3)

    if save_auc:
        out_path = os.path.join(model_dir, "avg_model_test_auc.dat")
        with open(out_path, 'w') as f:
            f.write(f"{avg_auc} += {stderr_auc}\n")
        print(f"saved AUC summary results to: {out_path}")

    if test_split_id is not None:
        print(f"Number of test sets containing this sequence: {found_in_test}")
    print("====================")
    print(f"Average test AUC for phase separation prediction: {avg_auc} ± {stderr_auc}")

    if test_split_id is not None and correct_preds:
        true_label = "Phase separates" if true_labels[0] == 1 else "Does not phase separate"
        avg_accuracy = round(mean(correct_preds), 2)
        stderr_accuracy = round(stdev(correct_preds) / sqrt(len(correct_preds)), 3) if len(correct_preds) > 1 else 0.0
        print(f"\nSequence test/train index = {test_split_id}")
        print(f"Ground truth label: {true_label}")
        print(f"Model accuracy for this sequence: {avg_accuracy} ± {stderr_accuracy}")
        print(f"True labels across test sets: {true_labels}")
        print(f"Predicted labels across test sets: {predicted_labels}")
    elif test_split_id is not None:
        print(f"\nSequence test/train index = {test_split_id} was never in any test set")

    print("====================")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model_results_dir", type=str, required=True, help="Directory containing model performance results")
    parser.add_argument("--seq_ps_file", type=str, required=True, help="Path to the file with sequence IDs and phase-separation labels")
    parser.add_argument("--split_sum", action='store_true', help="Whether the model is split-sum")
    parser.add_argument("--sequenceID", type=str, default=None, help="Optional sequence ID string to analyze, e.g., '0001' or '0001_1000'")
    parser.add_argument("--combo_name", type=str, default=None, help="If split-sum, specify the model combination name directory to use")
    parser.add_argument("--save_auc", action='store_true', help="If set, saves the average AUC and stderr to a .dat file")

    args = parser.parse_args()

    test_split_id = None
    if args.sequenceID is not None:
        test_split_id = get_test_split_id(args.seq_ps_file, args.sequenceID)

    if args.split_sum:
        if not args.combo_name:
            raise ValueError("--combo_name must be provided when --split_sum is used")
        analyze_split_sum(args.model_results_dir, args.combo_name, test_split_id, args.save_auc)
    else:
        analyze_non_split_sum(args.model_results_dir, test_split_id, args.save_auc)

if __name__ == "__main__":
    main()

