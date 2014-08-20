# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_Q_matrix_grid_mult_fq_3_carver
#PBS -j eo
#PBS -l nodes=3:ppn=8,walltime=00:30:00
#PBS -q debug
#PBS -A m1871

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data"
numProcs=24
maxl=10
beam_sig=0.064614317 # in radians.  HERA beam at 150 MHz
del_bl=23 # in nanoseconds.  HERA
sqGridSideLen=10
lowerFreq=50 # in MHz
upperFreq=90 # in MHz
deltaFreq=2 # in MHz


module load python/2.7.3
module load numpy
module load matplotlib
module load openmpi-gnu
module load mpi4py

# For running on my laptop or Carver
echo "Starting run now..."
date
mpirun -np $numProcs python-mpi "$codeLoc/mpi_Q_matrix_grid_mult_fq_3.py" $outputLoc $maxl $beam_sig $del_bl $sqGridSideLen $lowerFreq $upperFreq $deltaFreq
date