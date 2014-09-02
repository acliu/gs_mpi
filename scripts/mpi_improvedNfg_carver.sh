# !/bin/bash
#PBS -S /bin/bash
#PBS -N mpi_improvedNfg_carver
#PBS -j eo
#PBS -l nodes=1:ppn=8,walltime=00:30:00
#PBS -q debug
#PBS -A m1871

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi" #"/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi"
Gmatrix_codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/monte_carlo" # "/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/monte_carlo"
skyInfoDIR='/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates' #'/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/skyTemplates'
ClFname="Cl_empirical.dat"
templateFname="GSM_templateMap_70MHz_nside32.dat" #"template_small.dat"
GmatrixLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data_test" #"/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data_test" #"/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles"
KfgFname="$outputLoc/Kfg.npy"
nside=32 #8 #32
nlmax=200
lowerFreq=120 # in MHz
upperFreq=150 #122 # in MHz
deltaFreq=1 # in MHz
beam_sig=1.57 
del_bl=0.844 
sqGridSideLen=12 #2
variableBeam=0
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

numProcs=8

module swap PrgEnv-pgi PrgEnv-gnu
module load gcc
module load openmpi-gcc
module load python/2.7.3
module load numpy
module load matplotlib
module list

#gfortran generate_Kfg.f90 -o generate_Kfg.x -O3
#mpif90 generate_Kfg_mpi.f90 -o generate_Kfg_mpi.x -O3
#ftn -O3 -ffast-math -funroll-loops -o generate_Kfg_mpi.x generate_Kfg_mpi.f90
#mpif90 -O3 -ffast-math -funroll-loops -o mpi_generate_Kfg.x mpi_generate_Kfg.f90

pushd $skyInfoDIR
rm -f args.dat
touch args.dat
echo $nside >> args.dat
echo $nlmax >> args.dat
echo $ClFname >> args.dat
echo $templateFname >> args.dat
echo Kfg.raw >> args.dat

echo "Forming the foreground covariance matrix in image space..."
date
#./generate_Kfg.x args.dat
#mpirun -np 7 generate_Kfg_mpi.x args.dat
time mpirun -np $numProcs "$codeLoc/./mpi_generate_Kfg.x" args.dat
time python "$codeLoc/raw2npy.py" Kfg.raw $KfgFname
rm -f Kfg.raw
date
echo "...Done!"
popd


echo "Computing G (instrumental response) matrices..."
date
time python "$Gmatrix_codeLoc/G_matrix_grid_mult_fq_3.py" $nside $lowerFreq $upperFreq $deltaFreq $beam_sig $del_bl $sqGridSideLen $variableBeam $GmatrixLoc
date
echo "...Done!"

echo "Multiplying matrices together..."
date
time python "$codeLoc/sandwichGKG.py" $nside $KfgFname $GmatrixLoc $del_bl $sqGridSideLen $beam_sig $variableBeam $lowerFreq $upperFreq $deltaFreq $outputLoc
date
echo "Done!"


