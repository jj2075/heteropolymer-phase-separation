import os
import csv
import numpy as np
import argparse
import math

def parse_cutoff(value):
    """Parse cutoff distance as int or float, depending on input"""
    try:
        if '.' in value:
            return float(value)  # Return as float
        else:
            return int(value)  # Return as integer
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid cutoff distance: {value}")

def calculate_k(i, j):
    """calculate wavenumber |k| for a given mode (i, j)"""
    return math.sqrt(i**2 + j**2)

def process_b2_from_data_files(seq_ids, b2_value, seq_list_path, output_dir):
    output_path = os.path.join(output_dir, 'b2.npy')
    os.makedirs(output_dir, exist_ok=True)

    feature_data = []

    if b2_value is not None:
        # HP model: use same b2 value for all seq_ids
        feature_data = [(str(seq_id), float(b2_value)) for seq_id in seq_ids]
    else:
        # IDP model: parse b2 value from csv
        b2_map = {}
        with open(seq_list_path, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            for row in reader:
                seqid = str(row[0]).zfill(4)
                b2 = float(row[5])
                b2_map[seqid] = b2
        
        for seq_id in seq_ids:
            if seq_id not in b2_map:
                raise ValueError(f"b2 value not found for seq_id {seq_id}")
            feature_data.append((str(seq_id), b2_map[seq_id]))

    np.save(output_path, np.array(feature_data, dtype=[('seqid', 'U10'), ('value', 'f8')]))
    print(f"Saved b2 data to {output_path}")

def process_rg_mean_variance_from_data_files(seq_ids, data_dir, output_dir, feature_name, cutoff):
    """process rg or cmap mean/variance data and save as .npy"""
    output_path = os.path.join(output_dir, f'{feature_name}.npy')
    os.makedirs(output_dir, exist_ok=True)

    feature_data = []
    for seq_id in seq_ids:
        if feature_name == 'rg':
            file_path = os.path.join(data_dir, f'avg_{feature_name}-{seq_id}.dat')
        else:
            # Use cutoff for mean/variance
            file_path = os.path.join(data_dir, f'{feature_name}_contact_map_{cutoff}_{seq_id}.dat')

        if os.path.isfile(file_path):
            value = np.loadtxt(file_path)
            feature_data.append((seq_id, value))
        else:
            print(f"{feature_name} file not found for sequence id: {seq_id} with contact with cutoff {cutoff}")

    np.save(output_path, np.array(feature_data, dtype=[('seqid', 'U10'), ('value', 'f8')]))
    print(f"saved {feature_name} data to {output_path}")

def process_partial_sums_from_data_files(seq_ids, dat_dir, output_dir, nyquist_frequency, cutoff):
    """Read partial power sums from .dat files, compute 'above' sums, and save aggregated data as .npy arrays."""
    output_dir = os.path.join(output_dir, "partial_powers")
    os.makedirs(output_dir, exist_ok=True)

    # Dictionary to aggregate data for each (row, col)
    aggregated_data = {}

    for seq_id in seq_ids:
        dat_file = os.path.join(dat_dir, f"partial_power_contact_map_{cutoff}_{seq_id}.dat")
        if not os.path.isfile(dat_file):
            print(f"Partial power file not found for sequence id: {seq_id} with cutoff {cutoff}")
            continue

        try:
            # Skip the header row (row,col,cumulative_power)
            data = np.loadtxt(dat_file, delimiter=',', dtype=[('row', 'i4'), ('col', 'i4'), ('cumulative_power', 'f8')], skiprows=1)

            # Compute total cumulative power (last row and column at Nyquist frequency)
            total_power = data['cumulative_power'].max()

            data.sort(order=['row', 'col'])

            for row in range(nyquist_frequency + 1):
                for col in range(row, nyquist_frequency + 1):
                    # Find the cumulative power for the current (row, col)
                    current_power = data[(data['row'] == row) & (data['col'] == col)]['cumulative_power']
                    if len(current_power) == 0:
                        print(f"No power data found for (row={row}, col={col}) in sequence id: {seq_id}")
                        continue

                    current_power = current_power[0]

                    # Compute 'above' sum: total_power - current_power
                    above_sum = total_power - current_power

                    # Aggregate data for both sum types 'up_to_include' and 'above'
                    for power_type, partial_sum in {'up_to_include': current_power, 'above': above_sum}.items():
                        key = (row, col, power_type)
                        if key not in aggregated_data:
                            aggregated_data[key] = []
                        aggregated_data[key].append((seq_id, partial_sum))

        except ValueError as e:
            print(f"Error loading data from {dat_file}: {e}")

    # Save aggregated data for each (row, col, power_type)
    for (row, col, power_type), values in aggregated_data.items():
        filename = f"power_{power_type}_{row}_{col}.npy"
        filepath = os.path.join(output_dir, filename)
        np.save(filepath, np.array(values, dtype=[('seqid', 'U10'), ('partial_sum', 'f8')]))
        print(f"Saved aggregated partial power data ({power_type}) for (row={row}, col={col}) to {filepath}")

def main(seq_ids, b2_value, seq_list_path, rg_dir, mean_dir, variance_dir, partial_sums_dir, output_dir, nyquist_frequency, cutoff_distance):
    process_b2_from_data_files(seq_ids, b2_value, seq_list_path, output_dir)
    process_rg_mean_variance_from_data_files(seq_ids, rg_dir, output_dir, 'rg', cutoff=cutoff_distance)
    process_rg_mean_variance_from_data_files(seq_ids, mean_dir, output_dir, 'mean', cutoff=cutoff_distance)
    process_rg_mean_variance_from_data_files(seq_ids, variance_dir, output_dir, 'variance', cutoff=cutoff_distance)
    process_partial_sums_from_data_files(seq_ids, partial_sums_dir, output_dir, nyquist_frequency, cutoff=cutoff_distance)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="process features like b2, rg, variance, and partial sums")
    parser.add_argument("--seq_list_path", type=str, required=True, help="path to csv file containing sequence ids")
    parser.add_argument("--b2_value", type=int, default=None, help="global b2 value (for HP model only)")
    parser.add_argument("--rg_dir", type=str, required=True, help="directory containing rg data files")
    parser.add_argument("--mean_dir", type=str, required=True, help="directory containing contact map mean data files")
    parser.add_argument("--variance_dir", type=str, required=True, help="directory containing cmap variance data files")
    parser.add_argument("--partial_sums_dir", type=str, required=True, help="directory containing cumulative power files")
    parser.add_argument("--output_dir", type=str, required=True, help="directory to save npy array")
    parser.add_argument("--nyquist_frequency", type=int, default=10, help="set nyquist frequency to generate mode indices")
    parser.add_argument("--cutoff_distance", type=parse_cutoff, required=True, help="Cutoff distance for contact map (e.g., 2.5 or 16)")

    args = parser.parse_args()

    # read first col of CSV (skip header) 
    seq_ids = []
    with open(args.seq_list_path, 'r') as f:
        reader = csv.reader(f)
        next(reader)
        seq_ids = [str(row[0]).zfill(4) for row in reader if row] # pad seqID, e.g., 0004

    main(
        seq_ids,
        args.b2_value,
        args.seq_list_path,
        args.rg_dir,
        args.mean_dir,
        args.variance_dir,
        args.partial_sums_dir,
        args.output_dir,
        args.nyquist_frequency,
        args.cutoff_distance
    )
