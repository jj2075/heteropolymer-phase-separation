import argparse, math
import numpy as np
import scipy.integrate

def read_zbincounts(path):
    sections = {}
    N = 0
    with open(path, 'r') as f:
        for line in f:
            if len(line) > 1 and line[0] != '#':
                if N <= 0:
                    timestep = int(line.split()[0])
                    N = int(line.split()[1])
                    sections[timestep] = np.zeros(N)
                else:
                    sections[timestep][int(line.split()[0])-1] = float(line.split()[1])
                    N -= 1
    if N > 0:
        del sections[timestep]
    return sections

def center_profile(rho):
    zmax = np.argmax(rho)
    rho_centered1 = np.roll(rho, len(rho) // 2 - zmax)
    zcom = round(sum(i * rho_centered1[i] for i in range(len(rho))) / rho.sum())
    rho_centered2 = np.roll(rho_centered1, len(rho) // 2 - zcom)
    return rho_centered2

def calc_mean_centered_profile(sections):
    rho = np.zeros(len(sections[min(sections)]))
    rho2 = np.zeros(len(sections[min(sections)]))
    for profile in sections.values():
        rho_i = center_profile(profile)
        rho += rho_i
        rho2 += rho_i**2
    rho /= len(sections)
    std = np.sqrt(rho2 / len(sections) - rho**2)
    return rho, std

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str, help="path to log file")
    parser.add_argument('--neq', type=int, default=0, help="number of equilibration snapshots [0]")
    clargs = parser.parse_args()

    sections = read_zbincounts(clargs.path)
    mintimestep = sorted(sections)[clargs.neq]
    print("Timesteps: N = %d, min = %d, max = %d" % (len(sections), mintimestep, max(sections)))
    sections = {timestep : section for timestep,section in sections.items() if timestep >= mintimestep}
    mean, std = calc_mean_centered_profile(sections)
    print("Writing average_zbincounts.dat")
    with open('average_zbincounts.dat', 'w') as f:
        for i in range(len(mean)):
            f.write("%d %g %g\n" % (i, mean[i], std[i]))
