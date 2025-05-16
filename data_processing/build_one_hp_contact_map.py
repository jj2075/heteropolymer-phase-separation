import os
import numpy as np
import ast
import math, argparse
import scipy.optimize
from numpy import loadtxt
from scipy.integrate import cumtrapz
from scipy.interpolate import CubicSpline
from collections import defaultdict
import matplotlib.pyplot as plt

"""
Script for building a contact map from decorrelated two-chain snapshots using pairwise monomer distances. Snapshots are analyzed to identify monomer-monomer contacts based on a specified distance threshold, with contact probabilities accumulated across snapshots.
Writes to an output file the contact probabilities for each monomer pair, along with mean and standard deviation of contact probabilities.
The PMF fit function is used to filter snapshots based on the specific energy threshold before contributing to the contact map.
"""

def read_and_fit_pmf(file_path):
    # Read the pmf data from file
    data = np.loadtxt(file_path)
    r_values = data[:, 0]
    pmf_values = data[:, 1]
    # error_values = data[:, 2]
    
    # Fit a cubic spline to the data
    cs = CubicSpline(r_values, pmf_values, bc_type=('not-a-knot', 'clamped'))
    
    return cs

def read_timesteps(dump_path, com_path=None, gyr_path=None):
    # Read atom data from dump_path; or read atom data, com data, and gyr data from respective paths:
    def get_timestep_dump(f):
        aformat = 'iiifff' # (atomid, molid, type, x, y, z)
        data = {'box' : [], 'atoms' : {}}
        for line in f:
            if len(line) > 0 and line[0] != '#':
                if 'ITEM' in line:
                    item = line.split(':')[1][1:].strip()
                else:
                    if item == 'TIMESTEP':
                        timestep = int(line.strip())
                    elif item == 'NUMBER OF ATOMS':
                        data['natoms'] = int(line.strip())
                        data['coords'] = np.zeros((data['natoms'], 3))
                    elif 'BOX BOUNDS' in item:
                        data['box'].append(tuple(float(x) for x in line.split()))
                    elif 'ATOMS' in item:
                        l = tuple(int(x) if af == 'i' else float(x) \
                                  for x,af in zip(line.split(), aformat))
                        data['atoms'][l[0]] = (l[1], l[2])
                        data['coords'][l[0] - 1,:] = l[3:]
                        if len(data['atoms']) == data['natoms']:
                            break
        if not all(k in data for k in ['atoms', 'natoms', 'box']) \
           or len(data['atoms']) != data['natoms']:
            raise EOFError

        return timestep, data

    with open(dump_path, 'r') as dump_file:
        while True:
            try:
                timestep_coords, coords = get_timestep_dump(dump_file)
                yield timestep_coords, coords
            except EOFError:
                break

def get_contact_map(snapshot_dir, contact_distance, F_threshold):
    contact_map, norm = None, 0.
    minr, maxr = np.inf, 0.
    num_points = 1000
    Fthresh = min(min(pmf_fit_fcn(np.linspace(0.05, 16, num_points))) + F_threshold, 0.)

    print("Contact map F_threshold =", Fthresh)
 

    snapshot_files = [os.path.join(snapshot_dir, filename) for filename in os.listdir(snapshot_dir)]

    for path in snapshot_files:
        snapshots_path = []
        for timestep, coord_data in read_timesteps(path):

            mols = [np.array([coord_data['coords'][i,:] \
                                      for i in range(coord_data['coords'].shape[0]) \
                                  if coord_data['atoms'][i+1][0] == molid]) for molid in (1, 2)]
            coms = [np.mean(mol, axis=0) for mol in mols]
            dist = np.linalg.norm(coms[1] - coms[0])

            snapshots_path.append((dist, mols))

        snapshots = snapshots_path[:]

        if contact_map is None:
            contact_map = np.zeros((len(snapshots[0][1][0]), len(snapshots[0][1][1])))

        for dist, mols in snapshots:
            F = pmf_fit_fcn(dist)
            if F <= Fthresh:
                w = 1.
                if dist < minr:
                    minr = dist
                if dist > maxr:
                    maxr = dist
            else:
                w = 0.

            if w > 0.:
                norm += w
                for i in range(len(mols[0])):
                    for j in range(len(mols[1])):
                        if np.linalg.norm(mols[1][j, :] - mols[0][i, :]) < contact_distance:
                            contact_map[i, j] += w
                
    contact_map /= norm
    contact_map_mean = np.mean(contact_map)
    contact_map_std = np.std(contact_map)

    print("Contact map distance range: r = [%g, %g]" % (minr, maxr))
    print(f"Writing contact_map_{contact_distance}_{seq1id}_{seq2id}.dat")

    output_filename = f"contact_map_{contact_distance}_{seq1id}_{seq2id}.dat"
    output_path = os.path.join(clargs.output_folder, output_filename)

    with open(output_path, 'w') as f:
        f.write("# contact distance = %g\n" % contact_distance)
        f.write("# r = [%g, %g]\n" % (minr, maxr))
        f.write("\n# i j contact_probability mean_contact_probability stddev_contact_probability\n")
        for i in range(contact_map.shape[0]):
            for j in range(contact_map.shape[1]):
                f.write("%d %d %g %g %g\n" % (i, j, contact_map[i, j], contact_map_mean, contact_map_std))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('pmf_fit_path', type=str, help="path to PMF fit file")
    parser.add_argument('--bin_width', type=float, default=0.1, help="bin_width(sigma) in ABF simulation")
    parser.add_argument('--Rg_sum', type=float, default=2.51468, help="sum of the radius of gyrations of two chains far from each other ")
    parser.add_argument('--seq1id', type=str)
    parser.add_argument('--seq2id', type=str)
    parser.add_argument('--snapshot_dir', type=str, required=True, help="directory paths containing snapshot data")
    parser.add_argument('--contact_distance', type=float, default=3., help="threshold for determining a monomer-monomer contact [1.]")
    parser.add_argument('--F_threshold', type=float, default=1., help="threshold for defining PMF well [1.]")
    parser.add_argument('--output_folder', type=str, help="folder to write contact map to")
    clargs = parser.parse_args()

    filename = clargs.pmf_fit_path

    if not os.path.exists(filename):
      filename += ".BAK"

    bin_width = clargs.bin_width
    Rg_sum = clargs.Rg_sum
    seq1id = clargs.seq1id
    seq2id = clargs.seq2id

    # pmf_fit_fcn, pmf_fit, B2 = get_B2(filename, bin_width, Rg_sum)
    pmf_fit_fcn = read_and_fit_pmf(filename)

    # print("# B2 is {:.3f}".format(B2))
    print(f"Snapshot directory: {clargs.snapshot_dir}")

    if not os.path.isdir(clargs.snapshot_dir):
        raise FileNotFoundError(f"Snapshot directory '{clargs.snapshot_dir}' does not exist or is not a directory.")

    if clargs.snapshot_dir is not None:
        snapshot_dir = clargs.snapshot_dir
        print("snapshot_dir", snapshot_dir)
        contact_distance = clargs.contact_distance
        F_threshold = clargs.F_threshold

        contact_map = get_contact_map(snapshot_dir, contact_distance, F_threshold)
