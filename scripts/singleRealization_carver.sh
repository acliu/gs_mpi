# !/bin/bash

#PBS -S /bin/bash
#PBS -N singleRealization_carver
#PBS -j eo
#PBS -l nodes=1:ppn=1,walltime=00:30:00,pvmem=20GB
#PBS -q regular
#PBS -A m1871

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/monte_carlo"
#outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data/MCs"
outputLoc="/global/homes/a/acliu/globalSig/exactTemplateMatch/MCs"
#templateLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/haslam_408MHz_nside32.fits"
templateLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/gsm_haslam_32.fits"
#GmatrixLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data"
GmatrixLoc="/global/homes/a/acliu/globalSig/exactTemplateMatch"
pertFile="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/perturbationVariances.dat"

nside=32
beam_sig=0.31
del_bl=4.22
sqGridSideLen=12
lowerFreq=100 # in MHz
upperFreq=200 # in MHz
deltaFreq=1 # in MHz
numPertComponents=3
variableBeam=1
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz
echo "Starting run now..."
echo "Computing G matrices..."
date
python "$codeLoc/G_matrix_grid_mult_fq_3.py" $nside $lowerFreq $upperFreq $deltaFreq $beam_sig $del_bl $sqGridSideLen $variableBeam $GmatrixLoc
date
echo "...Done! Now creating directories for Monte Carlos"
date
python "$codeLoc/create_MC_directories.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $variableBeam
date
echo "...Done! Now doing the single realization"
python "$codeLoc/singleRealizationHaslamExtrap.py" $outputLoc $templateLoc $GmatrixLoc $nside $del_bl $sqGridSideLen $beam_sig $variableBeam $lowerFreq $upperFreq $deltaFreq $numPertComponents $pertFile
date
echo "...All Done!"

