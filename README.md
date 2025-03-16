## Thermodynamic integration error analysis paper 

This repository contains relevant scripts and notebooks for the research paper "Investigating errors in alchemical free energy predictions using random forest models and GaMD". 

`init_structures` contains the PDB structure for hu4d5-5 we used for all hu4d5-5 mutation predictions
`TI_scripts` contains scripts we used to generate Amber input files for thermodynamic integration.
`model_example` contains TI lambda production data for 13 hu4D5-5 cases, common functions used to process and plot the data, and a Jupyter notebook demonstrating our process in using GaMD end state sampling to mitigate the error due inadequate sampling of all 13 perturbations (if the model found a strong enough $R^2$ between the geometric and energetic DOF). The functions include

- `read_rotamers(runNum, rst=True)`
This function reads all the DV/DL values from the Amber TI output files for a directory of `runNum` runs, each containing 12 lambdas (weights and lambda values are determined by 12-point Gaussian quadrature). 

- Plotting functions for plotting 1-D and 2-D GaMD free energy profiles as well as their overlaying with TI lambda production data.

- Functions to read GaMD free energy profile PMF files output from PyReweighting-1D.py or PyReweighting-2D.py

- Jupyter notebooks contain the code and hyperparameters we used for our random forest model.


