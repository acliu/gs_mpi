# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_gsm_grid_mult_fq_3_carver
#PBS -j eo
#PBS -l nodes=3:ppn=8,walltime=00:15:00
#PBS -q regular
#PBS -A m1871

# Takes just 3 minutes and 14 seconds!

codeLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/"
outputLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles/"
templateMap='/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles/singleRealization/singleRealization_skyMaps_grid_del_bl_0.84_sqGridSideLen_2_fixedWidth_beam_sig_1.57.npy'
numProcs=4
maxl=2
nside=1
beam_sig=1.57 
del_bl=0.844 
sqGridSideLen=2
lowerFreq=120 # in MHz
upperFreq=122 # in MHz
deltaFreq=1 # in MHz
variableBeam=0
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

#beam_sig=0.064614317 # in radians.  HERA beam at 150 MHz
#del_bl=23 # in nanoseconds.  HERA


# For running on my laptop or Carver
echo "Starting run now..."
date
mpirun -np $numProcs python "$codeLoc/mpi_Nfg_singleRealization.py" $nside $outputLoc $lowerFreq $upperFreq $deltaFreq $maxl $beam_sig $del_bl $sqGridSideLen $variableBeam $templateMap
date

