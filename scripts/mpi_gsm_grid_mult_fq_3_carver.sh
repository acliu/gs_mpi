# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_gsm_grid_mult_fq_3_carver
#PBS -j eo
#PBS -l nodes=3:ppn=8,walltime=02:30:00
#PBS -q regular
#PBS -A m1871

# Takes just 3 minutes and 14 seconds!

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data"
templateMap='/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/GSM_templateMap_70MHz_nside32.fits'
numProcs=24
maxl=7
beam_sig=1.57 
del_bl=0.844 
sqGridSideLen=12
lowerFreq=120 # in MHz
upperFreq=150 # in MHz
deltaFreq=1 # in MHz

#beam_sig=0.064614317 # in radians.  HERA beam at 150 MHz
#del_bl=23 # in nanoseconds.  HERA

module load python/2.7.3
module load numpy
module load matplotlib
module load openmpi-gnu
module load mpi4py

# For running on my laptop or Carver
echo "Starting run now..."
date
mpirun -np $numProcs python-mpi "$codeLoc/mpi_gsm_grid_mult_fq_3.py" $templateMap $outputLoc $lowerFreq $upperFreq $deltaFreq $maxl $beam_sig $del_bl $sqGridSideLen
date

