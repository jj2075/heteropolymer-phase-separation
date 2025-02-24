# EOS simulations for IDPs

Example scripts and data files for running EOS simulations at fixed density to classify phase behavior.

## Files
- `start.lmp`: LAMMPS input script for equilibration
- `restart.lmp`: LAMMPS input script for production
- `sys.data`: LAMMPS data file
- `sys.settings`: LAMMPS force-field settings file

## Running simulations
1. **Equilibration**: Run `start.lmp` to perform equilibration. Estimate efficiency of simulation and determine simulation time required for `restart.lmp`
2. **Production**: Use `restart.lmp` for the longer production simulation; use data from the second half of the simulation for production pressures.

## Note
This model requires a potential not included in LAMMPS. For more details, please refer to the original paper:
Regy, R. M., Thompson, J., Kim, Y. C., & Mittal, J. (2021). Improved coarse‐grained model for studying sequence dependent phase separation of disordered proteins. Protein Science, 30(7), 1371-1379.
