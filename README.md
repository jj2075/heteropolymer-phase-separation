# Predicting heteropolymer phase separation using two-chain contact maps

Supporting data and scripts for:  
*Jessica Jin, Wesley Oliver, Michael A. Webb, William M. Jacobs (2025).*  
*Predicting Heteropolymer Phase Separation Using Two-Chain Contact Maps.*

## Data description

### HP model dataset
- **Sequences**: polymer sequences of length 20 with varying hydrophobicity patterns.
- **Phase-separation predictor data**: includes B_22, Rg, variance and mean of contact maps, and Sk values (power spectrum-derived statistics).
- **Contact maps**: Two-chain configurations averaged to generate 2D contact maps and related statistics.

### IDP model dataset
- **Global and challenge datasets**: Includes data for phase-separation predictors (B_22, Rg, variance and mean of contact maps, and Sk values). Includes contact maps for sequences in challenge dataset.

## Simulation and analysis scripts

- **HP model simulations**:
  - Run single-chain and two-chain simulations to calculate Rg and PMFs.
  - Generate contact maps and compute Fourier power spectra.

- **Logistic regression models**:
  - Train and evaluate logistic regression models using features derived from B22-matched datasets.
  - Analyze phase behavior predictors (*B22, Rg, contact map statistics*) for heteropolymer datasets.

- **Direct coexistence simulations**:
  - Create, equilibrate, and simulate slab geometries to analyze density profiles.
