import sys
import os
import numpy as n

beam_sigs = (n.pi/18,n.pi/6,5*n.pi/18,7*n.pi/18)
sqGridSideLens = (4,8,12,16)
variableBeams = (0,1)
lowerFreq = 100. # in MHz
lowerFreq /= 1000. # in GHz
numMCs = 10000
saveInterval = 150

kk=0
for beam_sig in beam_sigs:
    del_bl = 1/(2*n.pi*beam_sig*lowerFreq)
    print 'beam_sig = ',beam_sig,'; del_bl = ',del_bl
    for sqGridSideLen in sqGridSideLens:
        for variableBeam in variableBeams:
            fcontent = """# !/bin/bash

#PBS -S /bin/bash
#PBS -N MC_grid_{0}
#PBS -j eo
#PBS -l nodes=3:ppn=4,walltime=01:00:00,pvmem=5GB
#PBS -q regular
#PBS -A m1871

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/monte_carlo"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data/MCs"
templateLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/haslam_408MHz_nside32.fits"
GmatrixLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data"
pertFile="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/skyTemplates/perturbationVariances.dat"
nside=32
beam_sig={1}
del_bl={2}
sqGridSideLen={3}
lowerFreq=100 # in MHz
upperFreq=150 # in MHz
deltaFreq=1 # in MHz
numPertComponents=3
variableBeam={4}
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

numMCs={5}
saveInterval={6}
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
echo "...Done!"
echo "Consolidating MC results"
date
python "$codeLoc/consolidate_MCs.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $variableBeam $numMCs $saveInterval
date
echo "...all done!"
            
""".format(kk,beam_sig,del_bl,sqGridSideLen,variableBeam,numMCs,saveInterval)
            with open('./run_mpi_MC_grid_{0}.sh'.format(kk), 'w') as file:
                file.writelines(fcontent)
                # os.system('qsub run_mpi_MC_grid_{0}.sh'.format(kk))
            kk += 1

