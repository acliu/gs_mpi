# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_Q_matrix_grid_mult_fq_3_paper_carver
#PBS -j eo
#PBS -l nodes=9:ppn=8,walltime=04:00:00
#PBS -q regular
#PBS -A m1871

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data"

#codeLoc="/global/homes/m/mpresley/gs_mpi"
#outputLoc="/global/scratch2/sd/mpresley/gs_data"

numProcs=72
maxl=7 #30
del_bl=6.67 # for PAPER, assuming each antenna is ~2 meters wide #23 # in nanoseconds.  HERA
sqGridSideLen=5
lowerFreq=100 # in MHz
upperFreq=200 # in MHz
deltaFreq=1 # in MHz


module load python/2.7.3
module load numpy
module load matplotlib
module load openmpi-gnu
module load mpi4py

# For running on my laptop or Carver
echo "Starting run now..."
date
mpirun -np $numProcs python-mpi "$codeLoc/mpi_Q_matrix_grid_mult_fq_3_paper.py" $outputLoc $maxl $del_bl $sqGridSideLen $lowerFreq $upperFreq $deltaFreq
date
