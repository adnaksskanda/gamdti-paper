#!/bin/bash

in_parm=''
while getopts 'j:f:v' flag; do
  case ${flag} in
	j) jobname=${OPTARG} ;;
    f) in_parm=${OPTARG} ;;
    *) printf "Command not found. Enter -f for filename\n"
       exit 1
    esac
done

if [[ -z ${in_parm} ]]; then
  printf "Please use -f and provide a parm7 filename!\n"
  exit 1
fi

if [[ -z ${jobname} ]]; then 
  printf "Please use -j and provide a job name!\n"
  exit 1
fi

if [[ ! -f ./merge.out ]]; then
  printf "Please include merge.out file in the directory.\n"
  exit 1
fi

m_name=$(echo ${in_parm} | sed 's/......$//')


timask1=$(cat merge.out | grep "timask1")
timask2=$(cat merge.out | grep "timask2")
scmask1=$(cat merge.out | grep "scmask1")
scmask2=$(cat merge.out | grep "scmask2")



for i in $(seq 1 1 5); do 
#echo "$i"
mkdir run_${i}
cat > ./run_${i}/min.in << EOF
Minimization to relax initial bad contacts, explicit solvent
&cntrl
 imin=1, !Minimization
 maxcyc=5000, !maximum number of cycles
 ntmin=2,
 ntpr=1000, !print
 cut=8, !nonbond cut off
 ntr=1, !positional restraints
 restraintmask='@CA,C,N', !restraints on backbone CA, C, N
 restraint_wt=10.0, !weight for positional restraints
 
 iwrap=1,

 icfe=1, ifsc=1, clambda=0.5, scalpha=0.5, scbeta=12.0,
 ${timask1}
 ${timask2}
 ${scmask1}
 ${scmask2}
 clambda=0.5
&end
EOF


cat > ./run_${i}/heat.in << EOF
Explicit solvent heating
 &cntrl
  imin=0, !minimization off
  irest=0,ntx=1, !read coordinates, no velocities
  ntpr=1000, !energy information in mdout
  ntwx=1000, !coordinates will be written to trajectory
  nstlim=800000,!0.8ns heat run 0.15ns t=10 to 300K and 0.65ns t=300K nvt
  dt=0.001,
  ntt=3, !langevin dynamics
  gamma_ln=5.0, !collison frequency
  ig=-1,!random number generator
  ntc=1, 
  ntf=1, 
  cut=8, !nonbonded cutoff
  ntb=1, !constant volume
  tempi=10.0, !initial temperature
  temp0=300.0, !final temperature
  iwrap=1, !coordinates written to the restart and trajectory files will be "wrapped" into a primary box.
  ioutfm=1, !netcdf trajectory
  ntr=1, !positional restraint
  restraintmask='@CA,C,N',
  restraint_wt=1.0,

  icfe = 1, ifsc = 1,
  nmropt=1,
  ${timask1}
  ${timask2}
  ${scmask1}
  ${scmask2}
  clambda = 0.5,
 /
 &wt
  TYPE='TEMP0', ISTEP1=0, ISTEP2=150000,
  VALUE1=10.0, VALUE2=300.0,
 /
 &wt TYPE='END' /
EOF

cat > ./run_${i}/equil.in << EOF
Explicit solvent constant pressure equilibration
 &cntrl
  imin=0, !minimization off
  irest=1,ntx=5, !read coordinates,velocities
  ntpr=10000, !energy information in mdout
  ntwx=10000, !coordinates will be written to trajectory
  nstlim=4000000,!4ns equilibration run
  dt=0.001,
  ntt=3, !langevin dynamics
  gamma_ln=1.0, !collison frequency
  ig=-1,!random number generator
  ntc=1, 
  ntf=1, 
  cut=8, !nonbonded cutoff
  ntp=1, !constant pressure
  ntb=2, !when ntp=1
  barostat=2, !MC Barostat
  tempi=300.0, !initial temperature
  temp0=300.0, !final temperature
  iwrap=1, !coordinates written to the restart and trajectory files will be "wrapped" into a primary box.
  ioutfm=1, !netcdf trajectory
  ntr=1, !positional restraint
  restraintmask='@CA,C,N',
  restraint_wt=1.0,

  icfe = 1, ifsc = 1,
  ${timask1}
  ${timask2}
  ${scmask1}
  ${scmask2}
  clambda = 0.5,
 /
 &wt TYPE='END' /
EOF

cat > ./run_${i}/equilnorest.in << EOF
Explicit solvent constant volume equilibration
 &cntrl
  imin=0, !minimization off
  irest=1,ntx=5, !read coordinates,velocities
  ntpr=10000, !energy information in mdout
  ntwx=10000, !coordinates will be written to trajectory
  nstlim=5000000,!5 ns equilibration run
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
  clambda = 0.5,
 /
 &wt TYPE='END' /
EOF

cat > ./run_${i}/submission.sh << EOF
#!/bin/bash
#SBATCH -J ${jobname}_${i}
#SBATCH -p spot
#SBATCH --qos=spot
#SBATCH --gres=gpu:1

module purge
module load amber/20.12-intel-2021b-CUDA-11.4.1-AmberTools-20.15-Python-3.9.6

pmemd.cuda -O -i min.in -o min.out -c ../${m_name}.rst7 -p ../${m_name}.parm7 -r min.rst7 -ref ../${m_name}.rst7 -e min.en -inf min.mdinfo
pmemd.cuda -O -i heat.in -o heat.out -c min.rst7 -p ../${m_name}.parm7 -r heat.rst7 -ref min.rst7 -x TI_heating.nc -e heat.en -inf heat.mdinfo
pmemd.cuda -O -i equil.in -o equil.out -c heat.rst7 -p ../${m_name}.parm7 -r equil.rst7 -ref heat.rst7 -x TI_equil.nc -e equil.en -inf equil.mdinfo
pmemd.cuda -O -i equilnorest.in -o equilnorest.out -c equil.rst7 -p ../${m_name}.parm7 -r equilnorest.rst7 -ref equil.rst7 -x TI_equilnorest.nc -e equilnorest.en -inf equilnorest.mdinfo
EOF
done 
