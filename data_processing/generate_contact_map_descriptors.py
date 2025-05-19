import argparse
import os
import subprocess

"""
A driver script to compute and save contact map descriptors (mean, variance and Fourier outputs)
required for for phase separation prediction.

Example usage:
    For HP model:
        python3 generate_contact_map_descriptors.py --seqid 5 --model hp --b2_value 1000 --cutoff 3.0
    For IDP model:
        python3 generate_contact_map_descriptors.py --seqid 5 --model idp --cutoff 16
"""

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate statistical descriptors and Fourier outputs for phase separation prediction.")
    parser.add_argument('--seqid', type=int, required=True, help="Sequence ID (e.g., 5)")
    parser.add_argument('--model', type=str, required=True, choices=['hp', 'idp'], help="Model type: hp or idp")
    parser.add_argument('--b2_value', type=str, required=False, choices=['1000', '400'], help="B2 value for HP (required if model is HP)")
    parser.add_argument('--cutoff', type=float, required=True, help="Cutoff distance (float for HP, integer for IDP)")
    
    args = parser.parse_args()
  
    if args.model == 'hp' and not args.b2_value:
        parser.error("For HP model, the --b2_value arg is required (e.g., 1000 or 400)")
    return args

def get_cutoff_suffix_and_label(cutoff, model):
    if model == 'idp':
        cutoff_label = f"{int(round(cutoff))}"
        folder_suffix = f"cutoff_{cutoff_label}_angstroms"
    else:
        cutoff_label = f"{cutoff:.1f}"
        folder_suffix = f"cutoff_{cutoff_label}_sigma"
    return folder_suffix, cutoff_label

def get_output_paths(seqid, model, cutoff, b2_value=None):
    seqid_padded = f"{int(seqid):04d}"
    base_dir = "/home/jj2075/heteropolymer-phase-separation/heteropolymer-phase-separation"

    if model == 'hp':
        stat_model_subdir = f"hp_b2_{b2_value}"
    else:
        stat_model_subdir = "idp_challenge_set"

    folder_suffix, cutoff_label = get_cutoff_suffix_and_label(cutoff, model)

    contact_map_dir = os.path.join(base_dir, "contact_maps", stat_model_subdir, folder_suffix)
    contact_map_file = os.path.join(contact_map_dir, f"contact_map_{cutoff_label}_{seqid_padded}_{seqid_padded}.dat")
    mean_path = os.path.join(base_dir, "contact_map_statistics", "mean", stat_model_subdir, folder_suffix)
    var_path = os.path.join(base_dir, "contact_map_statistics", "variance", stat_model_subdir, folder_suffix)
    msd_path = os.path.join(base_dir, "contact_map_fourier_outputs", "explained_variances", stat_model_subdir, folder_suffix)
    partial_power_path = os.path.join(base_dir, "contact_map_fourier_outputs", "partial_power_sums", stat_model_subdir, folder_suffix)

    mean_file = os.path.join(mean_path, f"mean_contact_map_{cutoff_label}_{seqid_padded}.dat")
    var_file = os.path.join(var_path, f"variance_contact_map_{cutoff_label}_{seqid_padded}.dat")
    msd_file = os.path.join(msd_path, f"msd_contact_map_{cutoff_label}_{seqid_padded}.dat")
    power_sum_file = os.path.join(partial_power_path, f"partial_power_contact_map_{cutoff_label}_{seqid_padded}.dat")

    return contact_map_file, mean_file, var_file, msd_file, power_sum_file

def main():
    args = parse_arguments()

    contact_map_file, mean_file, var_file, msd_file, power_sum_file = get_output_paths(
        seqid=args.seqid, model=args.model, cutoff=args.cutoff, b2_value=args.b2_value
    )
    print("\n" + "=" * 80)  # separator
    # Call process_contact_map_statistics.py
    print("Processing contact map mean and variance...")
    subprocess.run([
        "python3", "process_contact_map_statistics.py",
        "--contact_map_path", contact_map_file, # path to the contact map
        "--mean_path", mean_file,               # path to save the mean of the contact map
        "--variance_path", var_file,            # path to save the variance of the contact map
        "--cutoff_distance", str(args.cutoff)   #cutoff distance
    ], check=True)

    print("-" * 80) # separator
    # Call process_contact_map_fourier_outputs.py
    print("Processing contact map FFT...")
    subprocess.run([
        "python3", "process_contact_map_fourier_outputs.py",
        "--contact_map_path", contact_map_file, # path to the contact map
        "--msd_output_path", msd_file,  # path to save the variance divergence index (MSD)
        "--cumulative_power_output_path", power_sum_file  # Path to save partial power sums
    ], check=True)

    print("-" * 80) # separator
    print(f"Successfully generated outputs for mean and variance:")
    print(f"  - Contact map mean file: {mean_file}")
    print(f"  - Contact map variance file: {var_file}")
    print("-" * 80) # separator
    print(f"Successfully generated Fourier-related outputs:")
    print(f"  - Variance divergence index (MSD) file: {msd_file} (not req'd for prediction)")
    print(f"  - Contact map partial power sums file: {power_sum_file}")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()
