# !/bin/bash

#PBS -S /bin/bash
#PBS -N logFiles/mpi_consolidate_MCs_carver
#PBS -j eo
#PBS -l nodes=1:ppn=1,walltime=00:20:00,pvmem=20GB
#PBS -q regular
#PBS -A m1871


codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/monte_carlo"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data/MCs"
beam_sig=1.57 
del_bl=0.844 
sqGridSideLen=12
lowerFreq=120 # in MHz
upperFreq=150 # in MHz
deltaFreq=1 # in MHz
variableBeam=0
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

numMCs=10000
saveInterval=150
numProcs=1

module load python/2.7.3
module load numpy
module load matplotlib

echo "Consolidating MC results"
date
python "$codeLoc/consolidate_MCs.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $variableBeam $numMCs $saveInterval
date
echo "...all done!"
