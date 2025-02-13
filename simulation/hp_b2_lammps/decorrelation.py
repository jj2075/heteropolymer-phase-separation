import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Calculate the center-of-mass (COM) distance decorrelation time')
parser.add_argument('seqid', type=str, help='Sequence ID for which to calculate decorrelation time')
args = parser.parse_args()

# file containing the last 500 lines of the colvars trajectory
filename = "last_500_traj.txt"

# Read COM distances
with open(filename, "r") as file:
    lines = [line for line in file.readlines() if not line.startswith("#")]
    com_distances = [float(line.split()[1]) for line in lines]

# Calculate the normed autocorrelation function
autocorr = np.correlate(com_distances, com_distances, mode='full')
autocorr = autocorr[autocorr.size // 2:]
autocorr /= autocorr[0]

time_step_increment = 10000  # colvars.traj output frequency (every 10,000 timesteps)
lag_times = np.arange(len(autocorr)) * time_step_increment

plt.plot(lag_times, autocorr, color='black')
plt.xlabel('Lag time')
plt.ylabel('Autocorrelation')
plt.title(f'Autocorrelation function for seqID {args.seqid}')

# Identify the decorrelation time (where autocorrelation falls below 1/e)
decorrelation_time_index = np.where(autocorr < 1 / np.e)[0][0]
decorrelation_time = decorrelation_time_index * time_step_increment
print(f"decorrelation time: {decorrelation_time} steps")
plt.axvline(x=decorrelation_time, color='red', linestyle='--', label='Decorrelation time')
plt.legend()
plt.close()

# Save the decorrelation time to be read later for sampling decorrelated configurations
with open('decorrelation_time.txt', 'w') as f:
    f.write(str(decorrelation_time))
