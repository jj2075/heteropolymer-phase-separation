#!/bin/bash
# Script to copy in-target decorrelated configurations based on the timesteps listed
# in winning_sequence.colvars.traj and store copied configurations in 'win-configs' directory.

mkdir -p win-configs

while read line; do
    # get the timestep (first column)
    timestep=$(echo "$line" | awk '{print int($1)}')

    # copy the corresponding configuration file to 'win-configs'
    config_file="configs/dump${timestep}.atom"
    if [[ -f "$config_file" ]]; then
        cp "$config_file" "win-configs/"
    else
        echo "Warning: $config_file not found. Skipping..."
    fi
done < winning_sequence.colvars.traj
