## Code and digital materials for Investigating errors in alchemical free energy predictions using random forest models and GaMD paper

This repository contains relevant scripts and notebooks for the research paper "Investigating errors in alchemical free energy predictions using random forest models and GaMD". 

`init_structures` contains the PDB structure for hu4d5-5 we used for all hu4d5-5 mutation predictions

`TI_scripts` contains scripts we used to generate Amber input files for thermodynamic integration.

`model_example` contains TI lambda production data for 13 hu4D5-5 cases, common functions used to process and plot the data, and a Jupyter notebook demonstrating our process in using GaMD end state sampling to mitigate the error due inadequate sampling of all 13 perturbations (if the model found a strong enough $R^2$ between the geometric and energetic DOF). Please note that the computed ddG value means and uncertainties may differ by a few hundredths digits from the published values due to the stochastic nature of bootstrapping. 

`model_example` contains subfolders `TI_data`, `common_functions`, `example_notebooks`, and `gamd_pmfs`
`TI_data` has the nearby rotamers, important interatomic distances, and energetic DV/DL values for each frame for our hu4D5-5 TI calculations.

`common_functions` provides utility functions to make the notebooks in `example_notebooks` readable. The utility functions provided in `common_functions` include:

- Functions to read all the DV/DL values from the Amber TI output files for a directory of `runNum` runs, each containing 12 lambdas (weights and lambda values are determined by 12-point Gaussian quadrature). 

- Plotting functions for plotting 1-D and 2-D GaMD free energy profiles as well as plotting these profiles overlaid with TI lambda production data.

- Functions to read GaMD free energy profile PMF files output from PyReweighting-1D.py or PyReweighting-2D.py

- Jupyter notebooks contain the code and hyperparameters we used for our random forest model.

`example_notebooks` contains a separate notebook for each TI ddG estimation, each notebook details our random forest + GaMD based TI correction approach for the specific mutation

`gamd_pmfs` contain the free energy profiles needed for the corrections in the `example_notebooks`

Please do not hesitate to reach out to me (listed as corresponding author in the manuscript) if you have any questions!


