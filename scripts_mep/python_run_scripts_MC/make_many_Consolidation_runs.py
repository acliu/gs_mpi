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
#PBS -N consolidate_grid_{0}
#PBS -j eo
#PBS -l nodes=1:ppn=1,walltime=00:20:00,pvmem=20GB
#PBS -q regular
#PBS -A m1871

codeLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_mpi/monte_carlo"
outputLoc="/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data/MCs"
beam_sig={1}
del_bl={2}
sqGridSideLen={3}
lowerFreq=100 # in MHz
upperFreq=150 # in MHz
deltaFreq=1 # in MHz
variableBeam={4}
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

numMCs={5}
saveInterval={6}
numProcs=1

module load python/2.7.3
module load numpy
module load matplotlib

echo "Consolidating MC results"
date
python "$codeLoc/consolidate_MCs.py" $outputLoc $del_bl $sqGridSideLen $beam_sig $lowerFreq $upperFreq $deltaFreq $variableBeam $numMCs $saveInterval
date
echo "...all done!"


""".format(kk,beam_sig,del_bl,sqGridSideLen,variableBeam,numMCs,saveInterval)
            with open('./run_mpi_consolidation_grid_{0}.sh'.format(kk), 'w') as file:
                file.writelines(fcontent)
            kk += 1

