#!/bin/bash

#author: Skanda Sastry
#purpose: To create 12-lambda folders for 5 replicates, each one having the right clambda flag, all other
#amber input flags identical. This is for thermodynamic integration after the structures have been fully
# minimized/heated/equilibrated.

# 12 lambdas for the 12-point Gaussian quadrature we will use to integrate
lambdas=(0.00922 0.04794 0.11505 0.20634 0.31608 0.43738 0.56262 0.68392 0.79366 0.88495 0.95206 0.99078)


if [[ ! -f ./merge.out ]]; then
  printf "Please include merge.out file in the directory.\n"
  exit 1
fi


timask1=$(cat merge.out | grep "timask1")
timask2=$(cat merge.out | grep "timask2")
scmask1=$(cat merge.out | grep "scmask1")
scmask2=$(cat merge.out | grep "scmask2")


for j in $(seq 1 1 5); do

cd run_${j}

for i in $(seq 0 1 11); do 
#echo "$i"
mkdir ${i}_lambda
cat > ${i}_lambda/pre_equil.in << EOF
NVT pre-prod equilibration
 &cntrl
  imin=0, !minimization off
  irest=0,ntx=1, !read coordinates,not velocities
  ntpr=10000, !energy information in mdout
  ntwx=10000, !coordinates will be written to trajectory
  nstlim=1000000,!1ns equilibration run
  dt=0.001,
  ntt=3, !langevin dynamics
  gamma_ln=1.0, !collison frequency
  ig=-1,!random number generator
  ntc=1, 
  ntf=1, 
  cut=8, !nonbonded cutoff
  ntb=1, !constant volume
  tempi=300.0, !initial temperature
  temp0=300.0, !final temperature
  iwrap=1, !coordinates written to the restart and trajectory files will be "wrapped" into a primary box.
  ioutfm=1, !netcdf trajectory

  icfe = 1, ifsc = 1,
  ${timask1}
  ${timask2}
  ${scmask1}
  ${scmask2}
  clambda = ${lambdas[${i}]},
 /
 &wt TYPE='END' /
EOF

cat > ${i}_lambda/lambda_prod.in << EOF
Explicit solvent constant volume production
 &cntrl
  imin=0, !minimization off
  irest=1,ntx=5, !read coordinates,velocities
  ntpr=25000, !energy information in mdout
  ntwx=25000, !coordinates will be written to trajectory
  ntwr=25000, !coordinates will be written to trajectory
  nstlim=5000000,!5ns production run
  dt=0.001,
  ntt=3, !langevin dynamics
  gamma_ln=2.0, !collison frequency
  ig=-1,!random number generator
  ntc=1, 
  ntf=1, 
  cut=8, !nonbonded cutoff
  ntb=1, !constant volume
  temp0=300.0, !final temperature
  iwrap=1, !coordinates written to the restart and trajectory files will be "wrapped" into a primary box.
  ioutfm=1, !netcdf trajectory
  
  icfe = 1, ifsc = 1,
  logdvdl = 1,
  ifmbar=1,
  mbar_states=14,
  mbar_lambda=0,0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738, 0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078,1,
  ${timask1}
  ${timask2}
  ${scmask1}
  ${scmask2}
  clambda = ${lambdas[${i}]},
 /
 &wt TYPE='END' /

EOF

done
cd ..
done
