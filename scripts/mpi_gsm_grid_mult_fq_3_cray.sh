# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_gsm_grid_mult_fq_3_cray
#PBS -j eo
#PBS -l mppwidth=24,walltime=00:30:00
#PBS -q debug
#PBS -A m1871

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data"
templateMap='/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/GSM_templateMap_70MHz_nside32.fits'
numProcs=24
maxl=10
beam_sig=0.064614317 # in radians.  HERA beam at 150 MHz
del_bl=23 # in nanoseconds.  HERA
sqGridSideLen=10
lowerFreq=120 # in MHz
upperFreq=160 # in MHz
deltaFreq=2 # in MHz

module load python/2.7.3
module load numpy
module load matplotlib
module load mpi4py

# For running on my laptop or Carver
echo "Starting run now..."
date
aprun -n 24 -N 24 python-mpi "$codeLoc/mpi_gsm_grid_mult_fq_3.py" $templateMap $outputLoc $lowerFreq $upperFreq $deltaFreq $maxl $beam_sig $del_bl $sqGridSideLen
date
