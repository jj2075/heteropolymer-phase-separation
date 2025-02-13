#!/usr/bin/env python3
# This script identifies the largest set of decorrelated center-of-mass (COM) configurations within a specified target range from a simulation trajectory
# It reads a .colvars.traj file containing timesteps and COM distances, then identifies the largest set of decorrelated configurations
# that fall within the target range. Targe configurations are saved to 'winning_sequence.colvars.traj'.

# The target range for COM distances can be adjusted by changing the `target_range_min` and `target_range_max` values below. The `decorrelation_time` is read from 'decorrelation_time.txt', which stores the minimum number of timesteps between decorrelated samples. This file was created from the previous LAMMPS simulation used to converge the PMF.

import argparse

parser = argparse.ArgumentParser(description='Process colvars.traj data to extract decorrelated COM-target configurations')
parser.add_argument('subdir_name', help='Name of the subdir containing the .colvars.traj file')
args = parser.parse_args()

subdir_name = args.subdir_name

# Load seq-specific decorrelation time
with open('decorrelation_time.txt', 'r') as file:
    decorrelation_time = float(file.read().strip())

# Define the target COM distance range
target_range_min = 1.0  # lower bound of the target COM distance range
target_range_max = 5.0  # upper bound of the target COM distance range

# Read timesteps and COM distances from the .colvars.traj file
timesteps = []
com_distances = []
with open(f'{subdir_name}.colvars.traj', 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue
        timestep, com_distance = map(float, line.split())
        timesteps.append(timestep)
        com_distances.append(com_distance)

winning_sequence = []  # list to store the largest sequence
winning_sequence_count = 0  # count of configs in the largest sequence

# Iterate over possible starting indices to find the largest set
for start in range(len(timesteps) // 3):  # limit starting index for efficiency
    last_timestep_within_range = -decorrelation_time
    sequence = []
    sequence_count = 0
    for i in range(start, len(timesteps)):
        timestep = timesteps[i]
        com_distance = com_distances[i]
        if target_range_min <= com_distance <= target_range_max:  # check if the COM distance is within the target range
            if timestep - last_timestep_within_range >= decorrelation_time:  # check decorrelation
                sequence.append((timestep, com_distance))
                sequence_count += 1
                last_timestep_within_range = timestep
    print(f"Number of decorrelated samples in COM dist range for sequence starting at index {start}: {sequence_count}")
    
    # Update the winning sequence (set) if a larger one is found
    if sequence_count > winning_sequence_count:
        winning_sequence_count = sequence_count
        winning_sequence = sequence

# Save the largest set of decorrelated COM-target configurations to 'winning_sequence.colvars.traj'
with open('winning_sequence.colvars.traj', 'w') as file:
    for timestep, com_distance in winning_sequence:
        file.write(f"{timestep:.1f} {com_distance:.3f}\n")

print(f"Largest set of decorrelated configs: {winning_sequence_count} samples")
