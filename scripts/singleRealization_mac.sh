# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_monte_carlo_gen_y_fq_3_carver
#PBS -j eo
#PBS -l nodes=3:ppn=4,walltime=00:30:00,pvmem=5GB
#PBS -q regular
#PBS -A m1871


codeLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/monte_carlo"
outputLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles/singleRealization"
templateLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/skyTemplates/gsm_haslam_1.fits"
GmatrixLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles/"
pertFile="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/skyTemplates/perturbationVariances.dat"
nside=1
beam_sig=1.57 
del_bl=0.844 
sqGridSideLen=2
lowerFreq=120 # in MHz
upperFreq=122 # in MHz
deltaFreq=1 # in MHz
numPertComponents=3
variableBeam=0
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

python "$codeLoc/create_MC_directories.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $variableBeam
date
python "$codeLoc/singleRealizationHaslamExtrap.py" $outputLoc $templateLoc $GmatrixLoc $nside $del_bl $sqGridSideLen $beam_sig $variableBeam $lowerFreq $upperFreq $deltaFreq $numPertComponents $pertFile $numMCs $saveInterval
date

