# !/bin/bash

codeLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/monte_carlo"
outputLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles/MCs"
templateLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/skyTemplates/haslam_408MHz_nside1.fits"
#DEBUG
GmatrixLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles"
pertFile="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/skyTemplates/perturbationVariances.dat"
nside=1 #DEBUG
beam_sig=0.064614317 # in radians.  HERA beam at 150 MHz
del_bl=23 # in nanoseconds.  HERA
sqGridSideLen=2
lowerFreq=50 # in MHz
upperFreq=90 # in MHz
deltaFreq=2 # in MHz
numPertComponents=3

numMCs=100
saveInterval=13
numProcs=4

python "$codeLoc/G_matrix_grid_mult_fq_3.py" $nside $lowerFreq $upperFreq $deltaFreq $beam_sig $del_bl $sqGridSideLen $GmatrixLoc

python "$codeLoc/create_MC_directories.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq

mpirun -np $numProcs python "$codeLoc/mpi_monte_carlo_gen_y_mult_fq_3.py" $outputLoc $templateLoc $GmatrixLoc $nside $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $numPertComponents $pertFile $numMCs $saveInterval

python "$codeLoc/consolidate_MCs.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $numMCs $saveInterval