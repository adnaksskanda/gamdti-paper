
This folder contains subfolders `TI_data`, `common_functions`, `example_notebooks`, and `gamd_pmfs`

`TI_data` has the nearby rotamers, important interatomic distances, and energetic DV/DL values for each frame for our hu4D5-5 TI calculations.

`common_functions` provides utility functions to make the notebooks in `example_notebooks` readable. The utility functions provided in `common_functions` include:

- Functions to read all the DV/DL values from the Amber TI output files for a directory of `runNum` runs, each containing 12 lambdas (weights and lambda values are determined by 12-point Gaussian quadrature). 

- Plotting functions for plotting 1-D and 2-D GaMD free energy profiles as well as plotting these profiles overlaid with TI lambda production data.

- Functions to read GaMD free energy profile PMF files output from PyReweighting-1D.py or PyReweighting-2D.py

- Functions to compute bootstrapping dGs given a dataframe with TI lambda production data 

- Code and hyperparameters we used for our random forest model are contained in `example_notebooks`

`example_notebooks` contains a separate notebook for each TI ddG estimation, each notebook details our random forest + GaMD based TI correction approach for the specific mutation

`gamd_pmfs` contain the free energy profiles needed for the corrections in the `example_notebooks`

`error_calculations` contains the data for the original TI predictions, the RF+GaMD corrections, and the extended 25 ns/lambda production predictions as well as python scripts which compute the mean and uncertainties for the RMSE and MAE for each of these sets of predictions.