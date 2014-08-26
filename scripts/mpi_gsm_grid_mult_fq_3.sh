# !/bin/bash

codeLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi"
outputLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles"
templateMap='/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/skyTemplates/GSM_templateMap_70MHz_nside2.fits'
numProcs=6
maxl=2
beam_sig=0.064614317 # in radians.  HERA beam at 150 MHz
del_bl=23 # in nanoseconds.  HERA
sqGridSideLen=2
lowerFreq=50 # in MHz
upperFreq=90 # in MHz
deltaFreq=2 # in MHz

# For running on my laptop or Carver
mpirun -np $numProcs python "$codeLoc/mpi_gsm_grid_mult_fq_3.py" $templateMap $outputLoc $lowerFreq $upperFreq $deltaFreq $maxl $beam_sig $del_bl $sqGridSideLen

