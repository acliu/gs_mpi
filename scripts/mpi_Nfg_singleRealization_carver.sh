# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_Nfg_singleRealization_carver
#PBS -j eo
#PBS -l nodes=3:ppn=8,walltime=00:15:00
#PBS -q regular
#PBS -A m1871

# Takes just 3 minutes and 14 seconds!

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi"
outputLoc="/global/homes/a/acliu/globalSig/exactTemplateMatch"
templateMap='/global/homes/a/acliu/globalSig/exactTemplateMatch/MCs/singleRealization_skyMaps_grid_del_bl_4.22_sqGridSideLen_12_lambdaBeam_beam_sig_0.31.npy'
numProcs=24
maxl=7
nside=32
beam_sig=0.31
del_bl=4.22
sqGridSideLen=12
lowerFreq=100 # in MHz
upperFreq=200 # in MHz
deltaFreq=1 # in MHz
variableBeam=1
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

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
mpirun -np $numProcs python "$codeLoc/mpi_Nfg_singleRealization.py" $nside $outputLoc $lowerFreq $upperFreq $deltaFreq $maxl $beam_sig $del_bl $sqGridSideLen $variableBeam $templateMap
date

