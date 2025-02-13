import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from scipy.interpolate import CubicSpline
from multiprocessing import Pool
import os
import time


def main():
    """
    Main function to generate and save contact maps for a set of sequences.

    Notes
    -----
    - Reads sequence indices from a CSV file.
    - Creates commands for each sequence and each cutoff distance.
    - Uses multiprocessing to calculate contact maps in parallel.
    - Saves the resulting contact maps to a NumPy file.
    """
    t0 = time.perf_counter()
    df  = pd.read_csv('/path/to/data_files/challenge_set.csv')
    all_indices = df['sequence_index'].values

    # Generate contact maps using these cutoff values 
    cutoffs = [8, 12, 16, 20, 24, 28, 32]

    print('generating commands')
    commands = []
    keys = []
    for idx in all_indices:
        for c in cutoffs:
            keys.append((idx, c))
            commands.append((idx, c, 0.4))
    print(f'{len(commands)} commands generated')

    print('calculating contact maps')
    # Use a Pool for parallel computation
    with Pool() as p:
        results = p.starmap(calc_contact_map, commands)
    print('contact maps calculated')

    print('saving contact maps')
    contact_maps = {}
    file_name = 'contact_maps.npy'
    # Associate each command key with its resulting contact map
    for key, contact_map in zip(keys, results):
        contact_maps[key] = contact_map
    np.save(file_name, contact_maps)
    print('contact maps saved')

    elapsed_time = seconds_to_hms(time.perf_counter()-t0)
    print(f'Time elapsed: {elapsed_time[0]:02}:{elapsed_time[1]:02}:{elapsed_time[2]:02}')

    return 0


def calc_contact_map(idx, cutoff, energy_offset=0.4):
    """
    Calculate a normalized contact map for a given sequence index and distance cutoff.

    Parameters
    ----------
    idx : int
        The sequence index for which to compute the contact map.
    cutoff : float
        The distance cutoff for considering two monomers in contact.
    energy_offset : float, optional
        An offset to shift the PMF (Potential of Mean Force) curve, by default 0.4kJ/mol.

    Returns
    -------
    numpy.ndarray
        A 2D array representing the average contact frequencey between
        monomers of sequences within the PMF well below 0.4kJ/mol.

    Notes
    -----
    - Navigates to subdirectories named "run..." within a parent directory
      related to the given sequence index.
    - Uses `get_pmf_zeros` to find zero points of the PMF, then leverages
      `pmf_com_correlation_time` to estimate the correlation time.
    - Accumulates contact maps from each subdirectory and divides
      by the total frame count to get a final normalized matrix.
    """
    parent = f'/path/to/simulation_scripts/poly{idx}'
    parent = f'/scratch/gpfs/wo6860/idp_b2/an_2024_seqs/b2_results/poly{idx}'
    subdirs = [os.path.join(parent, dir) for dir in os.listdir(parent) if dir.startswith('run')]

    # Retrieve the zero points of the PMF curve (two values)
    zeros = get_pmf_zeros(subdirs, energy_offset=energy_offset)

    # Initialize empty contact matrix
    sequence = seq_from_data(os.path.join(parent,'run0/sys.data'))
    contact_matrix = np.zeros((len(sequence), len(sequence)))

    # We'll keep track of the total number of frames that contributed
    count = 0

    for dir in subdirs:
        # Estimate correlation time from the COM distance time series
        tau = int(pmf_com_correlation_time(dir) + 1)
        temp_contact_matrix, temp_count = pmf_contact_map(dir, zeros[0], zeros[1], cutoff, step=tau)

        # Accumulate contact matrices and frame counts
        contact_matrix += temp_contact_matrix
        count += temp_count

    # Return the average contact map
    return contact_matrix / count


def pmf_contact_map(path, min_dist, max_dist, cutoff, step=1):
    """
    Calculate the contact map for a given protein structure and trajectory.

    Parameters
    ----------
    path : str
        The path to the directory containing `sys.data` and `coords.dcd`.
    min_dist : float
        The minimum distance threshold for counting frames (based on COM distance).
    max_dist : float
        The maximum distance threshold for counting frames (based on COM distance).
    cutoff : float
        The atom-atom distance cutoff used to mark a contact in the contact map.
    step : int, optional
        The stride (e.g., correlation time) for reading frames. Defaults to 1.

    Returns
    -------
    contact_map : numpy.ndarray
        A 2D array of shape (n_chain1, n_chain2) indicating how many times
        atoms i (in chain1) and j (in chain2) are within `cutoff`.
    count : int
        The number of frames (within [min_dist, max_dist]) that contributed
        to this contact map.

    Notes
    -----
    - Starts analysis at 25% of the trajectory for equilibration purposes.
    - For each included frame, calculates the center-of-mass (COM) distance
      between chain1 and chain2. If it lies within [min_dist, max_dist],
      we compute an atom-atom distance array and mark all distances < cutoff.
    - The final `contact_map` is a sum over those marked frames.
    """
    # Construct the file paths
    structure = path + '/sys.data'
    trajectory = path + '/coords.dcd'

    # Initialize the MDAnalysis universe
    u = mda.Universe(structure, trajectory,
                     atom_style='id resid type charge x y z vx vy vz',
                     topology_format='DATA',
                     format='DCD')

    # Get the box length from the universe
    L = u.dimensions[0]

    # We start from 25% into the trajectory to skip initial frames
    start = int(0.25 * len(u.trajectory))

    # Create arrays to store results
    delta_coms = np.zeros(len(u.trajectory[start::step]))
    count = 0

    # Select chain1 and chain2
    chain1 = u.select_atoms('resid 1')
    chain2 = u.select_atoms('resid 2')

    # Initialize an empty contact map of shape (len(chain1), len(chain2))
    contact_map = np.zeros((len(chain1), len(chain2)))

    # Loop over frames from 'start' onward, using a stride of 'step'
    for i, ts in enumerate(u.trajectory[start::step]):
        # Get positions relative to the first atom in each chain
        chain1_delta = chain1.positions - chain1[0].position
        chain2_delta = chain2.positions - chain2[0].position

        # Apply periodic boundary conditions (PBC)
        chain1_pbc = np.where(chain1_delta > L/2,
                              chain1_delta - L,
                              np.where(chain1_delta < -L/2, chain1_delta + L, chain1_delta))
        chain2_pbc = np.where(chain2_delta > L/2,
                              chain2_delta - L,
                              np.where(chain2_delta < -L/2, chain2_delta + L, chain2_delta))

        chain1_masses = chain1.masses
        chain2_masses = chain2.masses

        # Compute the center of mass (COM) for both chains
        chain1_com = np.sum(chain1_pbc.T * chain1_masses, axis=1) / np.sum(chain1_masses)
        chain2_com = np.sum(chain2_pbc.T * chain2_masses, axis=1) / np.sum(chain2_masses)

        # Calculate the COM distance
        delta_com = np.linalg.norm(chain1_com - chain2_com)
        delta_coms[i] = delta_com

        # Only consider frames where delta_com is in [min_dist, max_dist]
        if min_dist < delta_com < max_dist:
            count += 1
            # distance_array returns the pairwise distances between all
            # atoms in chain1 and chain2
            distance_matrix = distances.distance_array(chain1, chain2)
            # Mark positions where distance < cutoff
            temp = distance_matrix < cutoff
            # Convert booleans to int and accumulate
            contact_map += temp.astype(int)
    
    return contact_map, count


def pmf_com_correlation_time(path):
    """
    Compute an approximate correlation time from the center-of-mass distance
    between two chains in a trajectory.

    Parameters
    ----------
    path : str
        Path to the directory containing `sys.data` and `coords.dcd`.

    Returns
    -------
    float
        An estimate of the correlation time for the COM distance timeseries.

    Notes
    -----
    - Loads the MDAnalysis Universe from `sys.data` and `coords.dcd`.
    - Begins analysis at 25% into the trajectory to skip early frames.
    - Calls `correlation_time` on the resulting time series to get the final estimate.
    """
    # Construct the file paths
    structure = path + '/sys.data'
    trajectory = path + '/coords.dcd'

    # Initialize the universe
    u = mda.Universe(structure, trajectory,
                     atom_style='id resid type charge x y z vx vy vz',
                     topology_format='DATA', format='DCD')

    # Get the box length from the Universe
    L = u.dimensions[0]

    # We skip the first 25% of frames for equilibration
    start = int(0.25 * len(u.trajectory))

    # Prepare an array for storing the COM distance
    delta_coms = np.zeros(len(u.trajectory[start:]))

    # Define the two chains
    chain1 = u.select_atoms('resid 1')
    chain2 = u.select_atoms('resid 2')

    # Loop over frames
    for i, ts in enumerate(u.trajectory[start:]):
        chain1_delta = chain1.positions - chain1[0].position
        chain2_delta = chain2.positions - chain2[0].position

        # Apply PBC
        chain1_pbc = np.where(chain1_delta > L/2,
                              chain1_delta - L,
                              np.where(chain1_delta < -L/2, chain1_delta + L, chain1_delta))
        chain2_pbc = np.where(chain2_delta > L/2,
                              chain2_delta - L,
                              np.where(chain2_delta < -L/2, chain2_delta + L, chain2_delta))

        chain1_masses = chain1.masses
        chain2_masses = chain2.masses

        # COM of each chain
        chain1_com = np.sum(chain1_pbc.T * chain1_masses, axis=1) / np.sum(chain1_masses)
        chain2_com = np.sum(chain2_pbc.T * chain2_masses, axis=1) / np.sum(chain2_masses)

        delta_com = np.linalg.norm(chain1_com - chain2_com)
        delta_coms[i] = delta_com

    # Use correlation_time function to estimate correlation time
    return correlation_time(delta_coms)


def seconds_to_hms(seconds):
    """
    Convert a number of seconds into hours, minutes, and seconds.

    Parameters
    ----------
    seconds : float
        The total number of seconds to be converted.

    Returns
    -------
    tuple of int
        A tuple in the form (hours, minutes, seconds).

    Examples
    --------
    >>> seconds_to_hms(3661)
    (1, 1, 1)
    """
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    seconds = int(seconds % 60)
    return hours, minutes, seconds


def correlation_time(array):
    """
    Estimate the correlation time for a 1D time series using block averaging.

    Parameters
    ----------
    array : array_like
        A one-dimensional time series of data (e.g., distances).

    Returns
    -------
    float
        The correlation time (tau) estimated by block averaging.

    Notes
    -----
    - Splits the data into blocks of increasing size and checks how the
      block averages deviate from the overall mean.
    - Calls `block_tau` for each block size, then averages the final half
      of these estimates to get a stable result.
    """
    # Initialize values and arrays for calculation
    length = len(array)
    avg = np.average(array)
    base_err = np.std(array)
    block_lim = 20
    err_array = np.zeros(((int(length/block_lim))-1))

    # Calculate an estimate of correlation time for each block size
    for block in range(1, length//block_lim):
        err_array[block-1] = block_tau(array, length, block, avg, base_err)

    # Based on the calculated tau's, gets the limit of tau with increasing block size
    return np.average(err_array[-1*int(0.5*len(err_array)):])


def block_tau(array, length, block, avg, base_err):
    """
    Calculate an estimate of the correlation time using block averaging.

    Parameters
    ----------
    array : array_like
        The 1D time series data.
    length : int
        The total length of `array`.
    block : int
        The block size used for partial averaging.
    avg : float
        The overall mean of the data.
    base_err : float
        The base standard deviation of the underlying data.

    Returns
    -------
    float
        An estimated correlation time for the specified block size.

    Notes
    -----
    - Accumulates block averages over the entire array, compares them
      to the global mean, and sums the squares of the differences.
    - Normalizes the final value to produce a tau estimate.
    """
    tau = 0
    block_avg = 0

    # We'll only use a total of (length//block)*block points
    for i in range(int(length/block)*block):
        block_avg += array[i]
        if ((i+1) % block) == 0:
            block_avg /= block
            block_avg = block_avg - avg
            tau += block_avg * block_avg
            block_avg = 0

    tau /= (length//block)
    tau *= block
    tau /= base_err*2

    return tau


def seq_from_data(data_file):
    """
    Extract a polymer sequence from a LAMMPS data file.

    Parameters
    ----------
    data_file : str
        The path to the data file containing the structure (i.e., sys.data).

    Returns
    -------
    str
        The extracted sequence as a string of residue identifiers.

    Notes
    -----
    - Reads the number of atoms and atom types from the header lines.
    - Assumes a certain file layout with data starting at specific lines.
    - Maps integer atom types to residue labels, then reconstructs the sequence
      from the 'Atoms' section.
    """
    # Read in file
    with open(data_file, 'r') as f:
        lines = f.readlines()

    # Get relevant parameters from first few lines in data file
    n_atoms = int(lines[2].split()[0])
    atom_types = int(lines[3].split()[0])
    n_seqs = n_atoms - int(lines[4].split()[0])
    seq_length = n_atoms // n_seqs

    # Get the atoms as strings to access via index
    residues = [str(lines[i].split()[-1]) for i in range(13, atom_types+13)]

    # The "Atoms" section typically starts around line 16 + atom_types
    # We build the sequence by reading the 3rd column (type) for each line
    atoms_start_at = 16 + atom_types
    sequence = ''.join([residues[int(lines[i].split()[2])-1]
                        for i in range(atoms_start_at, atoms_start_at + seq_length)])

    return sequence


def get_pmf_zeros(subdirs, energy_offset=0.4):
    """
    Calculates the zeros of the PMF curve averaged over multiple independent
    ABF simulations.

    Parameters
    ----------
    subdirs : list of str
        List of subdirectories containing PMF data files named 'out.pmf'.
    energy_offset : float, optional
        Energy offset to be subtracted from the PMF values. Default is 0.4.

    Returns
    -------
    list of float
        A list of zeros of the fitted PMF curve within the range [0, 40].

    Notes
    -----
    - For each subdir, loads 'out.pmf' (r, u).
    - Takes the average across subdirs, then applies a 2*T*kb*log(r) correction.
    - Shifts the PMF values so they sit above a baseline plus the `energy_offset`.
    - Fits a cubic spline and finds the roots within [0, 40].
    """
    pmf_data = []
    for dir in subdirs:
        # Load the PMF data from the file 'out.pmf'
        pmf_data.append(np.loadtxt(os.path.join(dir, 'out.pmf')))
    pmf_data = np.array(pmf_data)
    pmf_data = np.mean(pmf_data, axis=0)
    # Extract r and u values from pmf_data
    r_values = pmf_data[:, 0]
    u_values = pmf_data[:, 1]

    T = 300
    kb = 0.0019872041

    # Apply 2*T*kb*log(r) correction
    u_values += 2*T*kb*np.log(r_values)

    # Zero the data assuming PMF goes to zero at large r
    u_values -= np.mean(u_values[-10:]) 
    u_values -= np.min(u_values) + energy_offset

    # Create the cubic spline object
    cs = CubicSpline(r_values, u_values)

    # If the spline at r=0 is negative, we add zero at r=0 plus a positive root
    if cs(0) < 0:
        zeros = [0]
        zeros.append(min(zero for zero in cs.roots() if zero > 0))
    else:
        zeros = cs.roots()

    # Limit the zeros to between 0 and 40
    zeros = [zero for zero in zeros if 0 <= zero <= 40]
    return zeros


if __name__ == '__main__':
    main()