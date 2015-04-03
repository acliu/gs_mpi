# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_consolidate_MCs_carver
#PBS -j eo
#PBS -l nodes=1:ppn=1,walltime=00:20:00,pvmem=20GB
#PBS -q regular
#PBS -A m1871


codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/monte_carlo"
outputLoc="/global/homes/a/acliu/globalSig/beamSizeVariations/fq_100_200_1.7radBeam/MCs"
beam_sig=1.7
del_bl=0.94
sqGridSideLen=12
lowerFreq=100 # in MHz
upperFreq=200 # in MHz
deltaFreq=1 # in MHz
variableBeam=1
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

numMCs=10000
saveInterval=20
numProcs=12

module load python/2.7.3
module load numpy
module load matplotlib

echo "Consolidating MC results"
date
python "$codeLoc/consolidate_MCs.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $variableBeam $numMCs $saveInterval
date
echo "...all done!"
