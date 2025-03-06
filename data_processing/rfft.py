import os
import numpy as np
import numpy.ma as ma
import argparse
from numpy.fft import ifft2
import matplotlib.pyplot as plt
from collections import OrderedDict

"""
    Script for analyzing symmetrized contact maps using the 2D discrete real Fourier transform (RFFT2).
    Routines include iteratively adding Fourier modes in order of increasing wavenumber to compute
    partial sums of the power spectrum, approximating the contact map via the inverse RFFT2,
    and computing the explained variance (EV) for each approximation, defined as S^k/S, where S^k
    is the cumulative power up to the kth unique mode, accounting for symmetries, and S is the total power.
    Additionally, the script computes the variance divergence index (VDI), defined as the mean-squared difference
    between the wavenumber-ordered EV and power-ordered EV distributions. The VDI is a measure of the deviation
    of a given power spectrum from a typical contact map power spectrum in which the lower wavenumber
    (longer wavelength) modes dominate. Script saves VDI values and cumulative powers for each sequence to output files. 
"""

def get_contact_map(contact_map_path):
    if os.path.exists(contact_map_path):
        data = np.loadtxt(contact_map_path, skiprows=4)
        i_values = data[:, 0].astype(int)
        j_values = data[:, 1].astype(int)
        contact_probabilities = data[:, 2]
        num_monomers = max(np.max(i_values), np.max(j_values)) + 1
        contact_map = np.zeros((num_monomers, num_monomers))
        for i, j, prob in zip(i_values, j_values, contact_probabilities):
            contact_map[i, j] = prob

        contact_map = (contact_map + contact_map.T ) / 2
        # Normalize so DC component doesn't dominate
        contact_map_normalized = contact_map - np.mean(contact_map)

    else:
        print("contact map path doesnt exist")
    return contact_map, contact_map_normalized

def get_contact_map_idp(contact_map_path):
    if os.path.exists(contact_map_path):
        contact_map = np.loadtxt(contact_map_path)

        contact_map = (contact_map + contact_map.T) / 2
        contact_map_normalized = contact_map - np.mean(contact_map)

    else:
        print("Contact map path doesn't exist")
        return None, None

    #print("contact_map_symmetrized=", contact_map)

    #print("contact_map_symmetrized_normed=", contact_map_normalized)
    return contact_map, contact_map_normalized

def fft2(contact_map_normalized):
    ''' compute 2d FFT for the mean-centered symmetrized contact map '''
    fourier_transform_n = np.fft.fft2(contact_map_normalized)
    power_spectrum_n = np.abs(fourier_transform_n) ** 2

    return fourier_transform_n, power_spectrum_n

def rfft2(contact_map_normalized):
    ''' compute 2d RFFT for the mean-centered symmetrized contact map '''
    fourier_transform_n = np.fft.rfft2(contact_map_normalized)
    power_spectrum_n = np.abs(fourier_transform_n) ** 2
    return fourier_transform_n, power_spectrum_n

def get_unique_rfft2_indices(N):
    """generate unique indices for a real, symmetric matrix (N x N).
    returns list of tuples"""
    unique_indices = [(0, 0)]  # include the DC component (0, 0)

    for col in range(1, N//2 + 1):
        rows = list(range(0, col + 1)) + list(range(N//2 + 1, N - col + 1))
        unique_indices.extend((row, col) for row in rows)

    return unique_indices

def get_unique_rfft2_indices_below_nyquist(N):
    ''' generate unique rfft2 indices for a diagonally symmetric matrix, where kx and ky are below the Nyquist frequency.
    These represent distinct frequency components in the lower triangle of the output array.'''
    unique_indices = [(0, 0)]  # include the DC component (0, 0)
    for col in range(1, N//2 + 1):
      rows = list(range(0, col + 1))
      unique_indices.extend((row, col) for row in rows)

    return unique_indices

def get_unique_rfft2_indices_above_nyquist(N):
    ''' generate unique rfft2 indices for a diagonally symmetric matrix, wehre kx is above the Nyquist frequency.
    These represent distinct frequency components in the upper triangle of the output array.'''
    unique_indices = []
    for col in range(1, N//2 + 1):
      rows = list(range(N//2 + 1, N - col + 1))
      unique_indices.extend((row, col) for row in rows)

    return unique_indices

def get_unique_rfft2_indices_and_frequency_mappings(N):
    ''' map redundant RFFT2 frequency indices to their corresponding physical counterparts.
    e.g.,  (19, 1) is the horizontally flipped version of (1, 1) in physical space. However,
    it represents a unique Fourier mode in the frequency domain due to how RFFT2 stores
    the real-valued signal's transformed coefficients. This mapping ensures that all distinct
    modes are indexed while associating them with physical wavenumbers (up to Nyquist frequency).
    '''
    unique_indices_below_nyquist = get_unique_rfft2_indices_below_nyquist(N)
    unique_indices_above_nyquist = get_unique_rfft2_indices_above_nyquist(N)

    unique_indices_and_mappings = {index: [] for index in unique_indices_below_nyquist}

    for index_above_nyquist in unique_indices_above_nyquist:
        row, col = index_above_nyquist

        if row > N // 2 and col <= N // 2:
            # this index corresponds to a horizontally flipped counterpart (col, N - row)
            counterpart_index = (col, N - row)
            unique_indices_and_mappings[counterpart_index].append(index_above_nyquist)
        elif row <= N // 2 and col > N // 2:
            print(f"warning: idx {index_above_nyquist} outside expected bounds")
        else:
            print(f"warning: idx {index_above_nyquist} outside expected bounds")

    return unique_indices_and_mappings

def k(kx, ky):
  """Compute the wavenumber k given kx and ky."""
  return np.sqrt(kx**2 + ky**2)

def order_dict_by_k(unique_indices_and_mappings):
    #print("unique_indices_and_mappings.items(): ", unique_indices_and_mappings.keys())
    """Order the dictionary by k."""
    sorted_items = sorted(unique_indices_and_mappings.items(), key=lambda x: k(x[0][0], x[0][1]))
    #print("unique_indices_and_mappings sorted: ", sorted_items)
    return OrderedDict(sorted_items)

def save_cumulative_powers(sorted_unique_indices, cumulative_powers, cumulative_power_output_path):
    """Save cumulative powers with corresponding sorted unique indices to a text file."""
    with open(cumulative_power_output_path, 'w') as f:
        f.write("row, col, cumulative power\n")
        for (row, col), power in zip(sorted_unique_indices, cumulative_powers):
            f.write(f"{row}, {col}, {power:.6e}\n")
    print(f"cumulative powers saved to {cumulative_power_output_path}")

def reconstruct_full_fft2_from_rfft2(full_fourier_transform, rfft2_matrix, original_matrix, sorted_unique_indices, msd_output_path):
    """Reconstruct the full FFT matrix from rfft2 output using Hermitian symmetry and diagonal symmetry.
    Also computes  beween relative explained variance for increasing wavenumber (k) and descending power on the same x-axis
    (mode index) to specific output path.
    """
    if isinstance(original_matrix, list):
        original_matrix = np.array(original_matrix)
    if isinstance(rfft2_matrix, list):
        rfft2_matrix = np.array(rfft2_matrix)

    contact_map_normed = original_matrix - np.mean(original_matrix)
    unique_indices_and_mappings = get_unique_rfft2_indices_and_frequency_mappings(len(rfft2_matrix))
    ordered_dict_unique_indices_and_mappings = order_dict_by_k(unique_indices_and_mappings)

    Nx, Ny = original_matrix.shape
    full_matrix = np.zeros((Nx, Ny), dtype=np.complex128)
    cumulative_fourier_matrix = np.zeros_like(full_matrix, dtype=np.complex128)

    cumulative_fourier_matrices = []  # For storing partial Fourier matrices
    summed_powers = [] # Tracking power values
    total_power = np.sum(np.abs(full_fourier_transform) ** 2) # i.e., S
    approximations = []  # Partial IFFT2 reconstructions
    mse_list = []  #  between partial reconstructions and original
    explained_variances_k_ordered = [] # explained variances ordered by wavenumber
    explained_variances_power_ordered = [] # explained variances ordered by power (descending)
    modes_used_for_partial_approximations = []  # Track modes used for each partial reconstruction
    modes_and_magnitudes = []  # To store each mode and its corresponding magnitude for sorting later

    def apply_symmetry(row, col, value):
        """Apply Hermitian, diagonal, composite, and boundary symmetry for a given index and return redundant indices."""
        redundant_indices = [(row, col)]
        full_matrix[row, col] = value
        cumulative_fourier_matrix[row, col] = value

        # General symmetry: Hermitian, diagonal, and composite
        if (row == 0 or col == 0) and not (row == 0 and col == 0):
            # Handle boundary symmetry
            if row == 0 and col != 0:
                # Apply diagonal symmetry
                full_matrix[col, row] = value
                cumulative_fourier_matrix[col, row] = value
                redundant_indices.append((col, row))

                # Apply Hermitian symmetry
                full_matrix[row, Ny - col] = np.conj(value)
                cumulative_fourier_matrix[row, Ny - col] = np.conj(value)
                redundant_indices.append((row, Ny - col))

                # Apply composite Hermitian + diagonal symmetry
                full_matrix[Ny - col, row] = np.conj(value)
                cumulative_fourier_matrix[Ny - col, row] = np.conj(value)
                redundant_indices.append((Ny - col, row))

            elif col == 0 and row != 0:
                # Apply diagonal symmetry
                full_matrix[col, row] = value
                cumulative_fourier_matrix[col, row] = value
                redundant_indices.append((col, row))

                # Apply Hermitian symmetry
                full_matrix[Nx - row, col] = np.conj(value)
                cumulative_fourier_matrix[Nx - row, col] = np.conj(value)
                redundant_indices.append((Nx - row, col))

                # Apply composite Hermitian + diagonal symmetry
                full_matrix[col, Nx - row] = np.conj(value)
                cumulative_fourier_matrix[col, Nx - row] = np.conj(value)
                redundant_indices.append((col, Nx - row))

        # Handle mapped physical counterparts (above Nyquist frequency)
        mapped_indices = ordered_dict_unique_indices_and_mappings.get((row, col), [])
        for mapped_row, mapped_col in mapped_indices:
            mapped_value = rfft2_matrix[mapped_row, mapped_col]

            full_matrix[mapped_row, mapped_col] = mapped_value
            cumulative_fourier_matrix[mapped_row, mapped_col] = mapped_value
            redundant_indices.append((mapped_row, mapped_col))

            # Apply Hermitian symmetry for mapped indices
            full_matrix[Nx - mapped_row, Ny - mapped_col] = np.conj(mapped_value)
            cumulative_fourier_matrix[Nx - mapped_row, Ny - mapped_col] = np.conj(mapped_value)
            redundant_indices.append((Nx - mapped_row, Ny - mapped_col))

            # Apply diagonal symmetry for mapped indices
            full_matrix[mapped_col, mapped_row] = mapped_value
            cumulative_fourier_matrix[mapped_col, mapped_row] = mapped_value
            redundant_indices.append((mapped_col, mapped_row))

            # Apply composite Hermitian + diagonal symmetry
            full_matrix[Ny - mapped_col, Nx - mapped_row] = np.conj(mapped_value)
            cumulative_fourier_matrix[Ny - mapped_col, Nx - mapped_row] = np.conj(mapped_value)
            redundant_indices.append((Ny - mapped_col, Nx - mapped_row))

        else:

            if row != 0 and col != 0:
                # Hermitian symmetry
                full_matrix[Nx - row, Ny - col] = np.conj(value)
                cumulative_fourier_matrix[Nx - row, Ny - col] = np.conj(value)
                redundant_indices.append((Nx - row, Ny - col))

                # Diagonal symmetry
                full_matrix[col, row] = value
                cumulative_fourier_matrix[col, row] = value
                redundant_indices.append((col, row))

                # Composite Hermitian + diagonal symmetry
                full_matrix[Ny - col, Nx - row] = np.conj(value)
                cumulative_fourier_matrix[Ny - col, Nx - row] = np.conj(value)
                redundant_indices.append((Ny - col, Nx - row))

                full_matrix[Nx - row, Ny - col] = np.conj(value)
                cumulative_fourier_matrix[Nx - row, Ny - col] = np.conj(value)
                redundant_indices.append((Nx - row, Ny - col))

            # Handle mapped physical counterparts (above Nyquist frequency)
            mapped_indices = ordered_dict_unique_indices_and_mappings.get((row, col), [])
            for mapped_row, mapped_col in mapped_indices:
                mapped_value = rfft2_matrix[mapped_row, mapped_col]

                full_matrix[mapped_row, mapped_col] = mapped_value
                cumulative_fourier_matrix[mapped_row, mapped_col] = mapped_value
                redundant_indices.append((mapped_row, mapped_col))

                # Apply Hermitian symmetry for mapped indices
                full_matrix[Nx - mapped_row, Ny - mapped_col] = np.conj(mapped_value)
                cumulative_fourier_matrix[Nx - mapped_row, Ny - mapped_col] = np.conj(mapped_value)
                redundant_indices.append((Nx - mapped_row, Ny - mapped_col))

                # Apply diagonal symmetry for mapped indices
                full_matrix[mapped_col, mapped_row] = mapped_value
                cumulative_fourier_matrix[mapped_col, mapped_row] = mapped_value
                redundant_indices.append((mapped_col, mapped_row))

                # Apply composite Hermitian + diagonal symmetry
                full_matrix[Ny - mapped_col, Nx - mapped_row] = np.conj(mapped_value)
                cumulative_fourier_matrix[Ny - mapped_col, Nx - mapped_row] = np.conj(mapped_value)
                redundant_indices.append((Ny - mapped_col, Nx - mapped_row))

        return redundant_indices

    for count, (row, col) in enumerate(ordered_dict_unique_indices_and_mappings.keys()):
        # Apply symmetry and collect all redundant indices for the current mode
        redundant_indices = apply_symmetry(row, col, rfft2_matrix[row, col])

        # Compute partial reconstruction using IFFT2
        partial_reconstruction_with_mean = np.real(ifft2(cumulative_fourier_matrix)) + np.mean(original_matrix)
        approximations.append(partial_reconstruction_with_mean)

        # Store the power of the mode for later sorting
        mode_magnitude = np.abs(rfft2_matrix[row, col]) ** 2
        modes_and_magnitudes.append(((row, col), mode_magnitude, redundant_indices))

        # Track modes used for this partial reconstruction
        modes_used_for_partial_approximation = [(row, col)] + redundant_indices
        modes_used_for_partial_approximations.append(modes_used_for_partial_approximation)

        # Sum power of unique modes in modes_used_for_partial_approximation
        #print("modes_used_for_partial_approximation ", modes_used_for_partial_approximation)
        unique_modes_for_partial_approximation = set(modes_used_for_partial_approximation)
        summed_power = np.sum([np.abs(cumulative_fourier_matrix[row, col]) ** 2 for row, col in unique_modes_for_partial_approximation])
        summed_powers.append(summed_power)

        # Compute S^k
        cumulative_power = np.sum(summed_powers)

        # compute EV = S^k / S
        explained_variance = (
          0.0 if (row == 0 and col == 0) else cumulative_power / total_power
        )
        explained_variances_k_ordered.append(explained_variance)

        # Compute mean squared error ()
        mse = np.mean((original_matrix - partial_reconstruction_with_mean) ** 2)
        mse_list.append(mse)

        if count == len(ordered_dict_unique_indices_and_mappings) - 1:
            print(f"Final mode {count}:  = {mse:.6f}")

        cumulative_fourier_matrices.append(cumulative_fourier_matrix.copy())

    # check that total cumulative power = total power of spectrum
    cumulative_powers = []
    #print("summed_powers: ", summed_powers)
    cumulative_sum = 0
    for power in summed_powers:
        cumulative_sum += power
        cumulative_powers.append(cumulative_sum)
    total_cumulative_power = cumulative_powers[-1]
    total_power = np.sum(np.abs(full_fourier_transform) ** 2)
    print(f"total_cumulative_power: {total_cumulative_power}, total_power: {total_power}")
    assert np.isclose(total_cumulative_power, total_power, atol=1e-6)

    # Final verification
    if not np.isclose(cumulative_fourier_matrix, full_fourier_transform, atol=1e-6).all():
        mismatch_indices = np.where(np.abs(cumulative_fourier_matrix - full_fourier_transform) > 1e-6)
        print(f"Mismatches at indices: {list(zip(mismatch_indices[0], mismatch_indices[1]))}")

    # Sort modes by magnitude (power) in descending order
    sorted_modes_by_power = sorted(modes_and_magnitudes, key=lambda x: x[1], reverse=True)

    # Compute explained variance for modes in descending order of power
    used_indices = set()

    for mode, magnitude, redundant_indices in sorted_modes_by_power:
        # only add if mode != (0,0)
        if mode == (0, 0):
            print("mode == (0,0), skipping")
            continue
        used_indices.add(mode)
        used_indices.update(redundant_indices)  # Add all redundant indices for this mode

        # Reconstruct the cumulative Fourier matrix using these indices
        temp_cumulative_fourier_matrix = np.zeros_like(full_matrix, dtype=np.complex128)
        for used_mode in used_indices:
            temp_cumulative_fourier_matrix[used_mode] = cumulative_fourier_matrix[used_mode]
        
        # compute S^k
        cumulative_power = np.sum(np.abs(temp_cumulative_fourier_matrix) ** 2)

        # compute EV
        explained_variance = cumulative_power / total_power
        explained_variances_power_ordered.append(explained_variance)

    # filter the (0,0) mode from relative_explained_variances
    explained_variances_k_ordered_filtered = explained_variances_k_ordered[1:]
    #print("len(explained_variances_k_ordered_filtered): ", len(explained_variances_k_ordered_filtered))
    #print("len(explained_variances_power_ordered): ", len(explained_variances_power_ordered))
    # Compute MSD between EV distributions (increasing k vs. descending power)
    msd_between_curves = np.mean((np.array(explained_variances_k_ordered_filtered) - np.array(explained_variances_power_ordered))**2)
    # print("MSD between EV(k-ordered) and EV(power-ordered): ", msd_between_curves)
    with open(msd_output_path, 'w') as msd_file:
        msd_file.write(f'{msd_between_curves:.6e}\n')

    return full_matrix, modes_used_for_partial_approximations, cumulative_fourier_matrices, cumulative_powers, explained_variances_k_ordered, explained_variances_k_ordered_filtered, explained_variances_power_ordered, msd_between_curves, approximations, mse_list

# ================================= PLOT EXPLAINED VARIANCE vs. k =================================
def plot_explained_variance(k_list, explained_variances_k_ordered):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    # Plot explained variance
    ax.plot(k_list, explained_variances_k_ordered, marker='o', linestyle='-', color='black', linewidth=1.7, markersize=6)
    ax.scatter(k_list, explained_variances_k_ordered, color='blue', edgecolor='#BEBEBE', s=60, zorder=3, marker='o', linewidths=0.8, clip_on=False)

    ax.set_xticks(np.arange(0, max(k_list) + 1, 2))
    ax.set_yticks(np.linspace(0, 1.0, 6))
    ax.set_ylim(0, 1.0)

    ax.set_xlabel(r'$k$', fontsize=14)
    ax.set_ylabel(r"EV $(\frac{\mathrm{Var}(X^k)}{\mathrm{Var}(X)})$", fontsize=14)

    ax.set_title("EV vs. $k$", fontsize=14, pad=10, loc="center")
    #ax.legend(loc="center right", fontsize=15, frameon=True)

    fig.subplots_adjust(top=0.8)

    # Style axis spines and ticks
    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_linewidth(1.0)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.xaxis.set_tick_params(labelsize=14)

    plt.tight_layout()
    #plt.show()

# ======================== PLOT EXPLAINED VARIANCE (increasing k and descending power) vs. mode index ==========================

def plot_explained_variance_comparison(explained_variances_k_ordered_filtered, explained_variances_power_ordered, msd_between_curves):
    """
    Plot explained variance for increasing wavenumber (k) and descending power on the same x-axis (mode index)
    """
    assert len(explained_variances_k_ordered_filtered) == len(explained_variances_power_ordered), "Lists must be same length"

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    indices = list(range(1, len(explained_variances_k_ordered_filtered) + 1))

    # Plot explained variance vs. increasing wavenumber (k)
    ax.plot(indices, explained_variances_k_ordered_filtered, marker='o', linestyle='-', color='blue', linewidth=1.7, markersize=6, label='$k$-ordered')
    ax.scatter(indices, explained_variances_k_ordered_filtered, color='blue', edgecolor='#B0C4DE', s=60, zorder=3, marker='o', linewidths=0.5)

    # Plot explained variance vs. descending power
    ax.plot(indices, explained_variances_power_ordered, marker='o', linestyle='-', color='green', linewidth=1.7, markersize=6, label='power-ordered (descending)')
    ax.scatter(indices, explained_variances_power_ordered, color='green', edgecolor='#98FB98', s=60, zorder=3, marker='o', linewidths=0.5)

    ax.set_xticks(np.arange(0, max(indices) + 1, 5))
    ax.set_yticks(np.linspace(0, 1.0, 6))
    ax.set_ylim(0, 1.0)

    ax.set_xlabel('Mode index (shared basis)', fontsize=14)
    ax.set_ylabel(r'Explained variance $(\frac{\mathrm{Var}(X^k)}{\mathrm{Var}(X)})$', fontsize=14)
    ax.set_title(f'Explained variance: $k$-ordered vs. power-ordered\nMSD: {msd_between_curves:.3e}', fontsize=16)

    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_linewidth(1.0)
    ax.yaxis.set_tick_params(labelsize=16)
    ax.xaxis.set_tick_params(labelsize=16)

    fig.subplots_adjust(top=0.8)

    ax.legend(fontsize=14)
    plt.tight_layout()
    #plt.show()

# ================================== PLOT PARTIAL FOURIER SPECTRA ==================================
def plot_partial_fourier_spectra(cumulative_fourier_matrices, full_fourier_transform, plot_indices):
    """Plot partial Fourier spectra (real parts of cumulative Fourier matrices) for modes used in approx."""
    fig, axes = plt.subplots(1, 4, figsize=(18, 3.5), gridspec_kw={'wspace': 0.2})
    titles = [r"FFT2$^{(1)}$", r"FFT2$^{(1)}$ + FFT2$^{(2)}$", r"FFT2$^{(1)}$ + $\dots$ + FFT2$^{(3)}$", "full FFT2"]

    all_real_parts = [cumulative_fourier_matrices[i].real for i in plot_indices]
    vmin, vmax = np.min(all_real_parts), np.max(all_real_parts)

    N = full_fourier_transform.shape[0]

    for i, ax in enumerate(axes):
        real_part = cumulative_fourier_matrices[plot_indices[i]].real

        colored_matrix = np.full((N, N), -1.0)
        np.putmask(colored_matrix, real_part != 0, real_part)
        masked_matrix = np.ma.masked_where(colored_matrix == -1, colored_matrix)

        im = ax.imshow(masked_matrix, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(titles[i], fontsize=14)

        ax.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
        ax.set_xticks(np.arange(0, N, 5))
        ax.set_xticklabels(np.arange(0, N, 5))
        ax.set_yticks(np.arange(0, N, 5))
        ax.set_yticklabels(np.arange(0, N, 5))

        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.75, pad=0.02, label="Fourier coefficient value")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    #plt.show()
    #plt.close()

# ================================= PLOT CMAP APPROXIMATIONS =================================
def plot_cmap_approximations(approximations, modes_used_for_partial_approximations, plot_indices):
    """Plot partial reconstructions (approximations) of the contact map."""
    fig, axes = plt.subplots(1, 4, figsize=(18, 3.5), gridspec_kw={'wspace': 0.2}, constrained_layout=True)
    titles = [r"$C^{(1)}$", r"$C^{(1)} + C^{(2)}$", r"$C^{(1)} + \dots + C^{(3)}$", "full matrix"]

    all_values = np.concatenate([approx.ravel() for approx in approximations])
    vmin, vmax = np.min(all_values), np.max(all_values)
    for i, ax in enumerate(axes):
        im = ax.imshow(approximations[plot_indices[i]], cmap='gray_r', vmin=vmin, vmax=vmax)
        ax.set_title(titles[i], fontsize=14)

        ax.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

        #ax.set_xticks([])
        #ax.set_yticks([])

        # print(f"Modes used for reconstruction: {modes_used_for_partial_approximations[plot_indices[i]]}")

    #plt.tight_layout()
    #plt.show()

# =========================================== DRIVER CODE ============================================
def main(seqid, contact_map_path, msd_output_path, cumulative_power_output_path):
    # Load contact map
    contact_map_symmetrized, contact_map_symmetrized_normed = get_contact_map(contact_map_path)

    # Compute FFTs and unique indices
    full_fourier_transform = fft2(contact_map_symmetrized_normed)[0]
    fourier_transform_n = rfft2(contact_map_symmetrized_normed)[0]
    unique_indices = get_unique_rfft2_indices(contact_map_symmetrized.shape[0])
    sorted_unique_indices = sorted(unique_indices, key=lambda x: np.sqrt(x[0]**2 + x[1]**2))

    # Reconstruct the Fourier transform and compute all outputs
    reconstructed_fft2, modes_used_for_partial_approximations, cumulative_fourier_matrices, cumulative_powers, \
    explained_variances_k_ordered, explained_variances_k_ordered_filtered, explained_variances_power_ordered, \
    msd_between_curves, approximations, mse_list = reconstruct_full_fft2_from_rfft2(
        full_fourier_transform=full_fourier_transform,
        rfft2_matrix=fourier_transform_n,
        original_matrix=contact_map_symmetrized,
        sorted_unique_indices=sorted_unique_indices,
        msd_output_path=msd_output_path
    )

    # Save cumulative powers
    unique_indices_below_nyquist = get_unique_rfft2_indices_below_nyquist(contact_map_symmetrized.shape[0])
    sorted_unique_indices_below_nyquist = sorted(unique_indices_below_nyquist, key=lambda x: np.sqrt(x[0]**2 + x[1]**2))
    #save_cumulative_powers(sorted_unique_indices_below_nyquist, cumulative_powers, cumulative_power_output_path)

    # Compute the full FFT2 for comparison
    difference = np.abs(full_fourier_transform - reconstructed_fft2)
    max_diff = np.max(difference)
    print(f"Max difference between full fft2 and reconstructed fft2 using unique modes: {max_diff:.5e}")

    # Plot explained variances and comparisons
    unique_indices_and_mappings = get_unique_rfft2_indices_and_frequency_mappings(len(fourier_transform_n))
    ordered_dict_unique_indices_and_mappings = order_dict_by_k(unique_indices_and_mappings)
    k_list = [k(x[0], x[1]) for x in ordered_dict_unique_indices_and_mappings.keys()]

    plot_explained_variance(k_list, explained_variances_k_ordered)
    plot_explained_variance_comparison(explained_variances_k_ordered_filtered, explained_variances_power_ordered, msd_between_curves)

    # Plot partial Fourier spectra and approximations
    plot_partial_fourier_spectra(cumulative_fourier_matrices, full_fourier_transform, plot_indices=[1, 2, 3, -1])
    plot_cmap_approximations(approximations, modes_used_for_partial_approximations, plot_indices=[1, 2, 3, -1])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process contact maps and compute FFT-related outputs")
    parser.add_argument('--seqid', type=str, required=True, help="seqID")
    parser.add_argument('--contact_map_path', type=str, required=True, help="Path to the contact map")
    parser.add_argument('--msd_output_path', type=str, required=True, help="Path to save the MSD between k-ordered EV and power-ordered EV")
    parser.add_argument('--cumulative_power_output_path', type=str, required=True, help="Path to save powers of Fourier approximations for seqID")

    args = parser.parse_args()

    main(
        seqid=args.seqid,
        contact_map_path=args.contact_map_path,
        msd_output_path=args.msd_output_path,
        cumulative_power_output_path=args.cumulative_power_output_path
    )
