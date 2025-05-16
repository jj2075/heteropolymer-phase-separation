import os
import argparse
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate mean and variance of symmetrized contact maps."
    )
    parser.add_argument("--contact_map_path", type=str, required=True, help="Path to the contact map")
    parser.add_argument("--mean_path", type=str, required=True, help="Path to save the contact map mean")
    parser.add_argument("--variance_path", type=str, required=True, help="Path to save the contact map variance")
    parser.add_argument("--cutoff_distance", type=parse_cutoff, required=True, help="Cutoff distance for contact map, e.g., 3.0 or 16")
    return parser.parse_args()

def parse_cutoff(value):
    """Parse cutoff distance as int or float, depending on input"""
    try:
        if '.' in value:
            return float(value)  # return as float
        else:
            return int(value)  # return as integer 
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid cutoff distance: {value}")

def get_contact_map(contact_map_path):
    """Load and compute a symmetrized contact map."""
    if not os.path.exists(contact_map_path):
        print("Contact map path doesn't exist:", contact_map_path)
        return None

    data = np.loadtxt(contact_map_path, skiprows=4)
    i_values, j_values, contact_probabilities = data[:, 0].astype(int), data[:, 1].astype(int), data[:, 2]
    num_monomers = max(np.max(i_values), np.max(j_values)) + 1
    contact_map = np.zeros((num_monomers, num_monomers))
    for i, j, prob in zip(i_values, j_values, contact_probabilities):
        contact_map[i, j] = prob
    return (contact_map + contact_map.T) / 2  # symmetrize

def save_statistic(statistic, output_path):
    """Save a single statistic (mean or variance) to a file"""
    # Ensure the directory exists before saving
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(f"{statistic}\n")

def main():
    args = parse_arguments()

    # Process contact map
    contact_map = get_contact_map(args.contact_map_path)
    if contact_map is not None:
        mean_value = np.mean(contact_map)
        variance_value = np.var(contact_map)

        # Save stats
        save_statistic(mean_value, args.mean_path)
        save_statistic(variance_value, args.variance_path)

        print(f"Mean saved to {args.mean_path}")
        print(f"Variance saved to {args.variance_path}")

if __name__ == "__main__":
    main()
