#!/bin/bash

# This script processes each seqID dir to read end of colvars.traj file and run decorrelation analysis.

sampling_configs_dir="<path_to_sampling_configs>"  # Dir containing sampling configuration scripts
data_dir="<path_to_data_files>"                   # Dir containing sequence and epsilon data
seqs_file="$data_dir/b2_400_combined_seqs-ps.lst" # File with sequence IDs

while read -r seq_id _; do
  cd "$seq_id" || continue

  traj_file=$(ls *.colvars.traj 2> /dev/null)
  if [[ -n "$traj_file" ]]; then
    tail -n 500 "$traj_file" > last_500_traj.txt
  fi

  python3 "$sampling_configs_dir/decorrelation.py" "$seq_id"

  cd ..
done < "$seqs_file"

