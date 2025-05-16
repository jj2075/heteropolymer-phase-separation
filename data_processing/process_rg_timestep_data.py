import os
import argparse
import pandas as pd

def compute_average_rg(rg_file_path):
    """
    Compute the average radius of gyration value (2nd column) from the given Rg file,
    excluding the first 10% of data 
    """
    data = []
    with open(rg_file_path, "r") as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.split()
                if len(parts) > 1:
                    data.append(float(parts[1]))

    total_lines = len(data)
    skip_count = int(total_lines * 0.1)
    data = data[skip_count:]

    return sum(data)/len(data)

def main():
    parser = argparse.ArgumentParser(description="compute average Rg from Rg-timestep data.")
    parser.add_argument("rg_dir", type=str, help="path to directory with Rg-timestep data")
    parser.add_argument("seqs_lst_path", type=str, help="path to sequences list file")
    parser.add_argument("save_dir", type=str, help="path to dir for saving avg_rg-{seqid}.dat files")
    args = parser.parse_args()

    # get sequence IDs
    seq_ids = []
    with open(args.seqs_lst_path, "r") as file:
        for line in file:
            if line.strip():
                seq_ids.append(line.split()[0]) # seqID in first col

    os.makedirs(args.save_dir, exist_ok=True)

    for seq_id in seq_ids:
        rg_file_path = os.path.join(args.rg_dir, seq_id, f"rg-{seq_id}.dat")
        save_file_path = os.path.join(args.save_dir, f"avg_rg-{seq_id}.dat")

        if not os.path.exists(rg_file_path):
            print(f"did not find {rg_file_path}, skipping.")
            continue

        avg_rg = compute_average_rg(rg_file_path)

        if avg_rg is not None:
            with open(save_file_path, "w") as save_file:
                save_file.write(f"{avg_rg}\n")
            print(f"saved mean Rg for {seq_id}") 
        else:
            print(f"No data data in {rg_file_path}")

if __name__ == "__main__":
    main()
