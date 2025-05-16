
# Predicting heteropolymer phase separation using two-chain contact maps

Supporting data and scripts for:  
**Jessica Jin, Wesley Oliver, Michael A. Webb, William M. Jacobs (2025).**  
*Predicting Heteropolymer Phase Separation Using Two-Chain Contact Maps.*

---

## Everything you need to predict phase separation from two-chain contact maps
Summary of key features of our end-to-end pipeline:
### 1. Generate contact-map features for logistic regression model training
For a given two-chain contact map, you can easily extract all necessary features for training logistic regression models.
See the [*data processing and feature extraction (`data_processing`)*](#data-processing-and-feature-extraction-data_processing) section below.

### 2. Train and evaluate logistic regression models
Test and train logistic regression models on polymer datasets.
See the [*logistic regression models (`logistic_regression`)*](#logistic-regression-models-logistic_regression) section below.

### 3. Evaluate accuracy of phase separation prediction by model
Assess how well logistic regression models models predict phase separation — at both the dataset level and per individual sequence. The pipeline supports:
- Evaluating test performance over randomized train/test splits for a polymer dataset
- Reporting a model's prediction accuracy for any given contact map
- Analyzing performance trends across models and preprocessing parameters (e.g., contact map cutoff distance)
See the [*logistic regression models (`logistic_regression`)*](#logistic-regression-models-logistic_regression) section for analysis tools.

---

## Data overview

### HP model dataset
- **Sequences**: length-20 heteropolymer sequences with varying hydrophobicity patterns
- **Features**: B_22, Rg, Var(C), Mean(C), S_k (Fourier-derived statistics)
- **Contact maps**: contact maps from two-chain simulations

### IDP dataset
- **`global_dataset.csv`**: 2034 sequences for the global IDP dataset with B_22 and Rg values
- **`idp_challenge_dataset.csv`**: Challenge dataset of 75 sequences with sequence-level features (e.g., charge fraction) and computed phase-separation predictors (B_22, Rg, Var(C), Mean(C), S_k), and contact maps
---

## simulation and analysis scripts

#### two-chain simulations for HP polymers (`simulation/hp_b2_lammps`)
- Run adaptive biasing force (ABF) simulations to calculate PMFs and compute B_22. The output `out.pmf` from the ABF simulation represents the bias applied during the ABF simulation; the true PMF is computed from this using `get_pmf.py`.

#### single-chain simulations for HP polymers (`simulation/hp_rg_lammps`)
- Compute Rg over time for each sequence to extract the average Rg used as a feature.

#### two-chain simulations for IDP polymers (`simulation/idp_b2_lammps`)
- Example scripts and data files for running two-chain simulations for IDPs to compute B2 and generate contact maps.

#### configuration sampling (`simulation/hp_sample_configurations_lammps`)
- Extract decorrelated configurations within a center-of-mass (COM) distance range for generating two-chain contact maps.

#### direct coexistence simulations (`simulation/hp_direct_coexistence_lammps`)
- Scripts for generating and analyzing slab geometries to classify phase behavior.

#### eos simulations (`simulation/idp_eos_lammps`)
- Example scripts and data files for running EOS simulations at fixed density to classify phase behavior.
---

## data processing and feature extraction (`data_processing`)
- **`build_hp_contact_map.py`** and **`build_idp_contact_map.py`**: Build two-chain contact maps from sampled configurations.  
- **`get_pmf.py`**: Calculate the PMF from `out.pmf` by removing the applied bias.  
- **`process_contact_map_statistics.py`**: Compute and save the mean and variance of contact maps.
- **`process_contact_map_fourier_outputs.py`**: Perform 2D Fourier decomposition of contact maps and compute and save the explained variance and partial sums of power spectrum.
- **`generate_contact_map_descriptors.py`**: The key driver script that takes a contact map as input and computes and saves contact map features/descriptors (mean, variance, explained variance, partial power sums) required for training a logistic regression model for phase separation prediction. This script calls `process_contact_map_statistics.py` and `process_contact_map_fourier_outputs.py` to perform computations.
- **`run_generate_contact_map_descriptors.sh`**: Bash script that automates the generation of contact map-based features for multiple contact maps by running `generate_contact_map_descriptors.py` for each contact map.
- **`consolidate_feature_arrays.py`**: Python script that aggregates feature data (e.g., mean, variance, partial sums) across sequences into feature-specific 2D NumPy arrays as {seq_id, feature_value}, streamlining data for model training.
- **`run_feature_array_consolidation.sh`**: a bash script that automates `consolidate_feature_arrays.py` for multiple datasets and features.
---

## logistic regression models (`logistic_regression`)

- **`model_feature_mappings/`**: Folder containing structured mappings that define which combinations of feature arrays (e.g., `rg.npy` and `b2.npy`) are used as inputs to specific logistic regression models (e.g., `"rg-b2"`).
- **`train_test_any_non_split_sum_model.py`**: Script for training and evaluating a logistic regression model using any feature combination (e.g., B2, Rg, contact map variance, etc.). Can be applied to both HP and IDP datasets using the specified feature set.
- **`train_test_a_split_sum_model.py`**: cript for training and evaluating a split-sum logistic regression model using any feature combination (e.g., B2, Rg, contact map variance, etc.). Applicable to any dataset where B2 is used as a feature, for either HP or IDP models; easily adaptable for B2-matched datasets.

### example job scripts for running models
- **`run_best_non_split_sum_model_hp_b2_mixed.sh`**:
Runs the best-performing non-split-sum logistic regression model on the HP B2-mixed dataset across 185 randomized train/test splits.
(This is not the best overall model reported in the paper. In the paper, results were reported for cutoff = 3.0 σ.)
- **`run_split_sum_model_hp_b2_mixed.sh`**:
Runs the split-sum logistic regression model on the HP B2-mixed dataset using power spectrum partitions across 185 randomized train/test splits.
(This IS the best-performing model for the HP B2-mixed dataset reported in the paper. In the paper, results were reported for cutoff = 3.0 σ.)
- **`run_best_non_split_sum_model_idp.sh`**:
Runs the best-performing non-split-sum model on the IDP challenge set (using four features: B2, single-chain Rg, contact-map mean, and contact-map variance).
(This performs equivalently to the split-sum model as the best-performing model for the IDP challenge set. In the paper, results were reported for cutoff = 24 Å.)

### analysis tools
- **`analyze_logistic_regression_model_results.py`**: Given a trained model, computes average classification accuracy across randomized test sets. Optional: for a given contact map ID, report  classification accuracy.
- **`run_analysis_of_logistic_regression_model_results.py`**: Wrapper script to automate model evaluation across multiple datasets and models.
- **`model_accuracy_vs_cutoff_distance.py`**: Evaluates how model performance changes as a function of the contact map cutoff distance—a useful sensitivity analysis.
---

## repository structure
- **`contact_maps/`**: Pre-computed contact maps for HP and IDP datasets at various cutoff contact distances. In this repo, HP contact maps include cutoff distances of 2.5σ, 3.0σ, and 3.5σ; IDP contact maps include distances of 20 Angstroms, 24 Angstroms, and 28 Angstroms. Other cutoff distances are excluded to save space but can be regenerated using the provided scripts.
- **`single_chain_rg/`**: Pre-computed single chain radius of gyration values for HP and IDP datasets. Values can be easily obtained by running `process_rg_timestep_data.py` in the `data_processing/` folder.
- **`phase_behavior_snapshots/`**: Displays snapshots of the final slab profiles from converged direct coexistence simulations for HP polymers.
  - These snapshots visualize the highest chain density above the chain overlap density at which phase behavior has converged. For exmaple, if a polymer exhibits phase separation at one chain density (computed as the number of chains divided by the total volume of the slab) but does not phase separate at a lower density that is above the chain overlap density, the polymer is classified as non-phase-separating, and the snapshot visualizes the highest density at which it does not exhibit phase separation.
  - Note: Phase behavior snapshots are not available for IDP sequences, which rely on density-based thermodynamic calculations for classification rather than direct coexistence simulations.
- **`sequences/`**: Sequence datasets and associated feature data.
  - Note: Sequence IDs between B2-matched datasets do not necessarily align. The datasets share most sequences; a few sequences that did not exhibit clear phase separation or aggregation in our direct-coexistence simulations at that B2 value were excluded from the dataset.
- **`simulation/`**: Job scripts and input files for running LAMMPS simulations.  
- **`data_processing/`**: Scripts for processing single-chain radius of gyration timestep data, computing contact maps and extracting Fourier-spectrum features.
- **`logistic_regression/`**: Logistic regression model training and evaluation.

---

## Reproducing results

1. **Run simulations**: use scripts in `simulation/` to run single-chain Rg and two-chain ABF simulations.  
2. **Process data**: use scripts in `data_processing/` to generate potential of mean force data, build contact maps from two-chain snapshot data, and compute contact-map features used for training logistic regression models.
3. **Train models**: use scripts in `logistic_regression/` to train and test logistic regression models.

---

## citation

If you use this repository, please cite:  
Jessica Jin, Wesley Oliver, Michael A. Webb, William M. Jacobs (2025).  
*Predicting Heteropolymer Phase Separation Using Two-Chain Contact Maps.*
