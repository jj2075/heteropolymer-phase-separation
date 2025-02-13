import argparse
import matplotlib.pyplot as plt
import numpy as np

# Load data from the density profile file
def load_density_data():
    filename = 'average_zbincounts.dat'
    try:
        return np.loadtxt(filename)
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: {filename} not found")

data = load_density_data()

bins = data[:, 0]
mean_counts = data[:, 1]
std_counts = data[:, 2]

# params for bin volume and density calculations
Lx, Ly, Lz = 20.0, 20.0, 180.0  # Box dimensions in sigma
Nbins = len(bins)
bin_volume = (Lx * Ly * Lz) / Nbins
length_of_polymer = 20

# Compute polymer densities and their standard deviations
polymer_densities = mean_counts / (bin_volume * length_of_polymer)
polymer_std_densities = std_counts / (bin_volume * length_of_polymer)

parser = argparse.ArgumentParser(description='Plot polymer chain density profile.')
parser.add_argument('--seqid', type=str, required=True, help='seqID (e.g., "0100")')
parser.add_argument('--epsilon', type=float, default=1.0, help='sequence epsilon value')
args = parser.parse_args()

# Plot the density profile
plt.errorbar(bins, polymer_densities, yerr=polymer_std_densities, fmt='-o', markersize=5)
title = f'Chain density profile for seqID {args.seqid} ($\epsilon$ = {args.epsilon}) using Lz = {Lz} $\sigma$'
plt.xlabel('z-bin index')
plt.ylabel('Polymer chain density ($\sigma^-3$)')
plt.ylim(0, 0.035)
plt.title(title)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'{args.seqid}-density-profile-Lz{int(Lz)}.pdf', dpi=300)
