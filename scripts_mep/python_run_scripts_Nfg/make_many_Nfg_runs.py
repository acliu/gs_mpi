import sys
import os
import numpy as n

beam_sigs = (n.pi/18,n.pi/6,5*n.pi/18,7*n.pi/18)
sqGridSideLens = (4,8,12,16)
variableBeams = (0,1)
lowerFreq = 100. # in MHz
lowerFreq /= 1000. # in GHz

kk=0
for beam_sig in beam_sigs:
    del_bl = 1/(2*n.pi*beam_sig*lowerFreq)
    print 'beam_sig = ',beam_sig,'; del_bl = ',del_bl
    for sqGridSideLen in sqGridSideLens:
        for variableBeam in variableBeams:
            fcontent = """# !/bin/bash

#PBS -S /bin/bash
#PBS -N Nfg_grid_{0}
#PBS -j eo
#PBS -l nodes=3:ppn=8,walltime=00:15:00
#PBS -q regular
#PBS -A m1871

# Takes just 3 minutes and 14 seconds!

codeLoc="/global/homes/m/mpresley/gs_mpi"
outputLoc="/global/scratch2/sd/mpresley/gs_data"
templateMap='/global/homes/m/mpresley/gs_mpi/skyTemplates/GSM_templateMap_70MHz_nside32.fits'

numProcs=24
maxl=7
beam_sig={1}
del_bl={2}
sqGridSideLen={3}
lowerFreq=100 # in MHz
upperFreq=150 # in MHz
deltaFreq=1 # in MHz
variableBeam={4}
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz

module load python/2.7.3
module load numpy
module load matplotlib
module load openmpi-gnu
module load mpi4py

# For running on my laptop or Carver
echo "Starting run now..."
date
mpirun -np $numProcs python-mpi "$codeLoc/mpi_gsm_grid_mult_fq_3.py" $templateMap $outputLoc $lowerFreq $upperFreq $deltaFreq $maxl $beam_sig $del_bl $sqGridSideLen $variableBeam
date
""".format(kk,beam_sig,del_bl,sqGridSideLen,variableBeam)
            with open('./run_mpi_Nfg_grid_{0}.sh'.format(kk), 'w') as file:
                file.writelines(fcontent)
                # os.system('qsub run_mpi_Nfg_grid_{0}.sh'.format(kk))
            kk += 1

