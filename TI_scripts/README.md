## Thermodynamic integration scripts

The shell scripts in these two folders are used to generate all the Amber input files necessary to do thermodynamic integration. Using a merged dual topology file (output from timerge in parmed), one can generate the input files using the commands listed below. Make sure that the merge.out file from tiMerge is included in the directory you execute these commands from.

``./struct_prep_replicates.sh -f {merged topology} -j {job name}``

``./lambda_maker.sh``
