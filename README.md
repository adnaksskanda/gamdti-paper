## Code and digital materials for Investigating errors in alchemical free energy predictions using random forest models and GaMD paper

This repository contains relevant scripts and notebooks for the research paper "Investigating errors in alchemical free energy predictions using random forest models and GaMD". 

`init_structures` contains the PDB structure for hu4D5-5 we used for all hu4D5-5 mutation predictions

`TI_scripts` contains scripts we used to generate Amber input files for thermodynamic integration.

`model_example` contains Jupyter notebooks demonstrating our random forest model trained on nearby side chain rotamers/interatomic distances and DV/DL. These are used to find the most influential degrees of freedom (DOF) for each mutation. The notebooks also demonstrate the corrections we made to our TI predictions by comparing them to the end state free energy profiles of the most influential DOF.  The folder also contains the necessary TI lambda production data to execute the model, the GaMD free energy profiles for the most influential DOFs for each mutation, and common utility functions used to process the data. Please note that the computed ddG value means and uncertainties may differ by a few hundredths digits from the published values due to the stochastic nature of bootstrapping. 


Please do not hesitate to reach out to me (listed as corresponding author in the manuscript) if you have any questions!

