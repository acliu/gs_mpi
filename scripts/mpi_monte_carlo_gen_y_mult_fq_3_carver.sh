# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_monte_carlo_gen_y_fq_3_carver
#PBS -j eo
#PBS -l nodes=1:ppn=8,walltime=00:30:00
#PBS -q debug
#PBS -A m1871


codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/monte_carlo"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data/MCs"
templateLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/haslam_408MHz_nside32.fits"
GmatrixLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data"
pertFile="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/perturbationVariances.dat"
nside=32
beam_sig=0.064614317 # in radians.  HERA beam at 150 MHz
del_bl=23 # in nanoseconds.  HERA
sqGridSideLen=5
lowerFreq=50 # in MHz
upperFreq=90 # in MHz
deltaFreq=2 # in MHz
numPertComponents=3

numMCs=100
saveInterval=11
numProcs=8

module load python/2.7.3
module load numpy
module load matplotlib
module load openmpi-gnu
module load mpi4py

echo "Starting run now..."
echo "Computing G matrices..."
date
python "$codeLoc/G_matrix_grid_mult_fq_3.py" $nside $lowerFreq $upperFreq $deltaFreq $beam_sig $del_bl $sqGridSideLen $GmatrixLoc
date
echo "...Done! Now creating directories for Monte Carlos"
date
python "$codeLoc/create_MC_directories.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq
date
echo "...Done! Now actually generating the Monte Carlos"
date
mpirun -np $numProcs python-mpi "$codeLoc/mpi_monte_carlo_gen_y_mult_fq_3.py" $outputLoc $templateLoc $GmatrixLoc $nside $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $numPertComponents $pertFile $numMCs $saveInterval
date
echo "...Done! Now consolidating the MC results"
date
python "$codeLoc/consolidate_MCs.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $numMCs $saveInterval
date
echo "...all done!"