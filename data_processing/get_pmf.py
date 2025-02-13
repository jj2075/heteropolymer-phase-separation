#/bin/python3
# calculate pmf and B2

import os
import numpy as np
import math, argparse
import scipy.optimize
from numpy import loadtxt
from scipy.integrate import cumtrapz
from scipy.interpolate import CubicSpline
from collections import defaultdict
import matplotlib.pyplot as plt

''' arguments: out_pmf_filename, bin_width, Rg_sum '''
''' Functions: correct and shift PMF, get B2 '''

def correct_data(filename, bin_width, Rg_sum):
  data = np.loadtxt(filename)

  ''' corrected data shape '''
  corrected_data = np.zeros((np.shape(data)[0]-1, np.shape(data)[1]))
  corrected_data[:,0] = data[:-1,0]+bin_width/2.0

  ''' Function: correction '''
  def gr_correction(r1, r2):
    if r1 > 0:
      return 2. * ((r2*np.log(r2) - r1*np.log(r1)) / (r2 - r1) - 1.)
    else:
      return 2. * (np.log(r2) - 1.)

  ''' Apply the correction '''
  dr = data[:,0]
  for i in range(np.shape(data)[0]-1):
      corrected_data[i,1] = data[i,1] + gr_correction(dr[i], dr[i+1])

  ''' Subtract the plateau '''
  def subtract_plateau(data_2d, Rg_sum):
    nearest_index = (np.abs(data_2d[:,0] - 2*Rg_sum)).argmin()
    avg_value = np.mean(data_2d[nearest_index:, 1])
    subtracted_data = data_2d[:, 1] - avg_value
    return np.column_stack((data_2d[:, 0], subtracted_data))

  ''' Return shifted, corrected data '''
  shifted_corrected_data = subtract_plateau(corrected_data, Rg_sum)
  # print("printing to pmf_corrected.txt... ")
  # np.savetxt(f'pmf_corrected-{seq1id}-{seq2id}.txt', shifted_corrected_data, delimiter=" ", fmt="%.6f")

  return shifted_corrected_data

def get_B2(filename, bin_width, Rg_sum):
  
  shifted_corrected_data = correct_data(filename, bin_width, Rg_sum)

  def mayer(f_array):
    ''' input: F_array (1D); returns: mayer_array (1D) = 1-exp(-F) '''
    y = np.zeros(np.shape(f_array))
    for i in range(np.shape(f_array)[0]):
      y[i] = 1.0 - np.exp(-f_array[i])
    return y

  def r2_mayer(r_array,f_array):
    y = np.zeros(np.shape(f_array))
    for i in range(np.shape(f_array)[0]):
      y[i] = r_array[i]**2 * (1. - np.exp(-f_array[i]))
    return y

  def pi_r2_mayer_function(r_array,f_array):
    y = np.zeros(np.shape(f_array))
    for i in range(np.shape(f_array)[0]):
      y[i] = 2.0 * math.pi * r_array[i]**2 * (1. - np.exp(-f_array[i]))
    return y

  def B2_from_mayer(r, b2_bins):
    return cumtrapz(b2_bins, r)

  def write_to_file(x_array, y_array, filename):
      with open(filename,'w') as file:
          for i in range(len(x_array)):
              file.write("{} {}\n".format(x_array[i], y_array[i]))

  ''' Calculate Mayer model '''
  data = np.loadtxt(filename)
  pmf_err = np.zeros((data.shape[0]-1, data.shape[1]))
  pmf_err[:,0] = data[:-1,0]+bin_width/2.0

  ''' compute the variance in a single PMF by using the variance of the data after 2Rg_sum '''
  def pmf_error(shifted_corrected_pmf_2d, Rg_sum):
    r_values = shifted_corrected_pmf_2d[:,0]
    # Find the index of the first value of dr_mean after 2 times Rg_sum
    idx = next((i for i, x in enumerate(r_values) if x > 2 * Rg_sum), None)

    for i in range(shifted_corrected_pmf_2d.shape[0]-1):
      if r_values[i] > 2 * Rg_sum:
        pmf_err[i, 1] = 24 * np.abs(shifted_corrected_pmf_2d[i, 1])
  
    return pmf_err[:,1]

  ''' compute B2 variance given a single r value '''
  def pmf_variance(i, shifted_corrected_pmf_2d, pmf_error_values):
      # pmf is the shifted
      r_values = shifted_corrected_pmf_2d[:,0]
      y_values = shifted_corrected_pmf_2d[:,1]
      return (r_values[i]**2 * np.exp(-y_values[i]) * pmf_error_values[i])**2

  def prior_variance(i, corrected_pmf_2d, Rg_sum):
      maxv = max(np.fabs(corrected_pmf_2d[:,1]))
      ''' The Gaussian envelope is centered around the maximum value (amplitude) maxv, and its width is controlled by gyrsuminf, which is the interaction range,
      which is approximated as the sum of the radii of gyration of the two interacting chains.
      The resulting value is multiplied by maxv to give the height of the Gaussian envelope.
      The expression represents a Gaussian distribution with a mean of 0 and a standard deviation of gyrsuminf. '''
      r_value = corrected_pmf_2d[:,0][i]
      gaussian_envelope_r = maxv * np.exp(-1. * (r_value / Rg_sum)**2)
      return (r_value **2 * (np.exp(gaussian_envelope_r) - 1.))**2
  
  ''' PMF error is variance past 2 * Rg_sum '''
  pmf_error_values = pmf_error(shifted_corrected_data, Rg_sum)

  ''' variance of data '''
  var_data = np.array([pmf_variance(i, shifted_corrected_data, pmf_error_values) for i in range(len(shifted_corrected_data[:,0]))])

  ''' variance model '''
  var_model = np.array([prior_variance(i, shifted_corrected_data, Rg_sum) for i in range(len(shifted_corrected_data[:,0]))])

  alpha = 0.1
  mayer_model = r2_mayer(shifted_corrected_data[:,0], shifted_corrected_data[:,1]) / (1. + alpha * (var_data / var_model))

  ''' get PMF cubic spline '''
  def get_f(r_value, mayer_model_value):
    return -math.log(-mayer_model_value / r_value**2 + 1.)

  pmf_model = np.array([get_f(shifted_corrected_data[:,0][i], mayer_model[i]) for i in range(len(shifted_corrected_data[:,0]))])

  for i in range(len(shifted_corrected_data[:,0])):
    if math.fabs(pmf_model[i]) < 0.04:
          pmf_model[i] = 0.

  ''' calculate mayer_data from pmf_model (which has smoothed out noise) '''
  mayer_data_smooth = pi_r2_mayer_function(shifted_corrected_data[:,0], pmf_model)
  B2_value = B2_from_mayer(shifted_corrected_data[:,0], mayer_data_smooth)[-1]

  ''' write mayer data and B2 data '''
  print("printing to mayer_data.txt...")
  mayer_filename = 'mayer_data-{}-{}.txt'.format(seq1id, seq2id)
  write_to_file(shifted_corrected_data[:,0], mayer(shifted_corrected_data[:,1]), mayer_filename)

  print("printing to B2_bins.txt...")
  bins_filename = 'B2_bins-{}-{}.txt'.format(seq1id, seq2id)
  write_to_file(shifted_corrected_data[:,0], mayer_data_smooth, bins_filename)

  ''' fit cubic spline for plotting '''
  pmf_model_cs = CubicSpline(shifted_corrected_data[:,0], pmf_model, bc_type=('not-a-knot', 'clamped'))

  def cs_extrap(r):
    if r < shifted_corrected_data[:,0][-1]:
        return pmf_model_cs(r)
    else:
        return 0.

  npoints_fine = 300
  def predict_fine(cs):
    return (np.array([r for r in np.linspace(0., shifted_corrected_data[:,0][-1], npoints_fine)]), \
      np.array([cs(r) for r in np.linspace(0., shifted_corrected_data[:,0][-1], npoints_fine)]))

  pmf_fit_fine = predict_fine(pmf_model_cs)

  ''' write PMF fit '''
  print("printing to pmf-fine.txt...")

  with open(f'pmf-fine-{seq1id}-{seq2id}.txt', 'w') as f:
       f.write("# B2 = %g\n" % B2_value)
       for i in range(len(pmf_fit_fine[0])):
                f.write("%g %g\n" % (pmf_fit_fine[0][i], pmf_fit_fine[1][i]))

  return B2_value

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('path', type=str, help="path to out.pmf file")
  parser.add_argument('--bin_width', type=float, default=0.1, help="bin_width(sigma) in ABF simulation")
  parser.add_argument('--Rg_sum', type=float, default=2.51468, help="sum of the radius of gyrations of two chains far from each other (Angstrom) ")
  parser.add_argument('--seq1id', type=str)
  parser.add_argument('--seq2id', type=str)

  clargs = parser.parse_args()
  filename = clargs.path
  bin_width = clargs.bin_width
  Rg_sum = clargs.Rg_sum
  seq1id = clargs.seq1id
  seq2id = clargs.seq2id

  B2 = get_B2(filename, bin_width, Rg_sum)
  print("# B2 is {:.3f}".format(B2))

# TO RUN: python3 get_pmf.py path --bin_width --Rg_sum
