# !/bin/bash

#PBS -S /bin/bash
#PBS -N logFiles/mpi_monte_carlo_gen_y_fq_3_carver
#PBS -j eo
#PBS -l nodes=3:ppn=4,walltime=00:30:00,pvmem=5GB
#PBS -q regular
#PBS -A m1871


codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/monte_carlo"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data/MCs"
templateLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/haslam_408MHz_nside32.fits"
GmatrixLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data"
pertFile="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/perturbationVariances.dat"
nside=32
beam_sig=1.57 
del_bl=0.844 
sqGridSideLen=12
lowerFreq=120 # in MHz
upperFreq=150 # in MHz
deltaFreq=1 # in MHz
numPertComponents=3
variableBeam=0
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

#beam_sig=0.064614317 # in radians.  HERA beam at 150 MHz
#del_bl=23 # in nanoseconds.  HERA
#sqGridSideLen=5
#lowerFreq=50 # in MHz
#upperFreq=90 # in MHz
#deltaFreq=2 # in MHz
#numPertComponents=3

numMCs=10000
saveInterval=150
numProcs=12

module load python/2.7.3
module load numpy
module load matplotlib
module load openmpi-gnu
module load mpi4py

echo "Starting run now..."
echo "Computing G matrices..."
date
python "$codeLoc/G_matrix_grid_mult_fq_3.py" $nside $lowerFreq $upperFreq $deltaFreq $beam_sig $del_bl $sqGridSideLen $variableBeam $GmatrixLoc
date
echo "...Done! Now creating directories for Monte Carlos"
date
python "$codeLoc/create_MC_directories.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $variableBeam
date
echo "...Done! Now actually generating the Monte Carlos"
date
mpirun -np $numProcs python-mpi "$codeLoc/mpi_monte_carlo_gen_y_mult_fq_3.py" $outputLoc $templateLoc $GmatrixLoc $nside $del_bl $sqGridSideLen $beam_sig $variableBeam $lowerFreq $upperFreq $deltaFreq $numPertComponents $pertFile $numMCs $saveInterval
date
echo "...all done!"
