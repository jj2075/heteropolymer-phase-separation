import os
import numpy as np
import argparse

def calculate_k(i, j):
    return np.sqrt(i**2 + j**2)

def generate_modes(nyquist_frequency):
    """generate unique modes up to the nyquist frequency"""
    modes = [(i, j) for i in range(nyquist_frequency + 1) for j in range(i, nyquist_frequency + 1)]
    return sorted(modes, key=lambda mode: calculate_k(*mode))

def combine_partial_sums(source_dir_400, source_dir_1000, output_dir, nyquist_frequency):
    """process partial sum data from two dirs"""
    partial_sums_output_dir = os.path.join(output_dir, "partial_sums")
    os.makedirs(partial_sums_output_dir, exist_ok=True)

    # generate modes up to the nyquist frequency
    modes_to_save = generate_modes(nyquist_frequency)

    for i, j in modes_to_save:
        for prefix in ["up_to_include", "above"]:
            file_name = f"power_{prefix}_{i}_{j}.npy"
            file_path_400 = os.path.join(source_dir_400, file_name)
            file_path_1000 = os.path.join(source_dir_1000, file_name)

            if os.path.exists(file_path_400) and os.path.exists(file_path_1000):
                data1 = np.load(file_path_400, allow_pickle=True)
                data2 = np.load(file_path_1000, allow_pickle=True)

                # add the b2 value as a suffix to the seqID
                data1_with_suffix = np.array([(seq + "_400", value) for seq, value in data1], dtype=data1.dtype)
                data2_with_suffix = np.array([(seq + "_1000", value) for seq, value in data2], dtype=data2.dtype)

                # combine the data
                combined_data = np.concatenate((data1_with_suffix, data2_with_suffix), axis=0)

                # save the combined data
                combined_file_path = os.path.join(partial_sums_output_dir, file_name)
                np.save(combined_file_path, combined_data)
                print(f"combined partial sums saved to {combined_file_path}")

def combine_simple_features(source_dir_400, source_dir_1000, output_dir, feature_name):
    """combine simple feature arrays like b2, rg, variance"""
    feature_output_dir = os.path.join(output_dir, feature_name)
    os.makedirs(feature_output_dir, exist_ok=True)

    file_name = f"{feature_name}.npy"
    file_path_400 = os.path.join(source_dir_400, file_name)
    file_path_1000 = os.path.join(source_dir_1000, file_name)

    if os.path.exists(file_path_400) and os.path.exists(file_path_1000):
        data1 = np.load(file_path_400, allow_pickle=True)
        data2 = np.load(file_path_1000, allow_pickle=True)

        # add the b2 value as a suffix to the seqID
        data1_with_suffix = np.array([(seq + "_400", value) for seq, value in data1], dtype=data1.dtype)
        data2_with_suffix = np.array([(seq + "_1000", value) for seq, value in data2], dtype=data2.dtype)

        # combine the data
        combined_data = np.concatenate((data1_with_suffix, data2_with_suffix), axis=0)

        # save the combined data
        combined_file_path = os.path.join(feature_output_dir, file_name)
        np.save(combined_file_path, combined_data)
        print(f"combined {feature_name} saved to {combined_file_path}")

def main(b2_dir_400, b2_dir_1000, rg_dir_400, rg_dir_1000, mean_dir_400, mean_dir_1000, variance_dir_400, variance_dir_1000, partial_sums_dir_400, partial_sums_dir_1000, output_dir, nyquist_frequency):
    combine_simple_features(b2_dir_400, b2_dir_1000, output_dir, "b2")
    combine_simple_features(rg_dir_400, rg_dir_1000, output_dir, "rg")
    combine_simple_features(mean_dir_400, mean_dir_1000, output_dir, "mean")
    combine_simple_features(variance_dir_400, variance_dir_1000, output_dir, "variance")
    combine_partial_sums(partial_sums_dir_400, partial_sums_dir_1000, output_dir, nyquist_frequency)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="combine feature arrays from two b2 values")
    parser.add_argument("--b2_dir_400", type=str, required=True, help="path to b2 directory for b2-400")
    parser.add_argument("--b2_dir_1000", type=str, required=True, help="path to b2 directory for b2-1000")
    parser.add_argument("--rg_dir_400", type=str, required=True, help="path to rg directory for rg-400")
    parser.add_argument("--rg_dir_1000", type=str, required=True, help="path to rg directory for rg-1000")
    parser.add_argument("--mean_dir_400", type=str, required=True, help="path to cmap mean directory for mean-400")
    parser.add_argument("--mean_dir_1000", type=str, required=True, help="path to cmap mean directory for mean-1000")
    parser.add_argument("--variance_dir_400", type=str, required=True, help="path to variance directory for variance-400")
    parser.add_argument("--variance_dir_1000", type=str, required=True, help="path to variance directory for variance-1000")
    parser.add_argument("--partial_sums_dir_400", type=str, required=True, help="path to partial power sums directory for b2-400")
    parser.add_argument("--partial_sums_dir_1000", type=str, required=True, help="path to partial power sums directory for b2-1000")
    parser.add_argument("--output_dir", type=str, required=True, help="directory to save combined data")
    parser.add_argument("--nyquist_frequency", type=int, default=10, help="set nyquist frequency for mode generation")

    args = parser.parse_args()
    main(
        args.b2_dir_400,
        args.b2_dir_1000,
        args.rg_dir_400,
        args.rg_dir_1000,
        args.mean_dir_400,
        args.mean_dir_1000,
        args.variance_dir_400,
        args.variance_dir_1000,
        args.partial_sums_dir_400,
        args.partial_sums_dir_1000,
        args.output_dir,
        args.nyquist_frequency
    )
