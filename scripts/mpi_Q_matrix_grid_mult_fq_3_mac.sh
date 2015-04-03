# !/bin/bash

#PBS -S /bin/bash
#PBS -N mpi_Q_matrix_grid_mult_fq_3_carver
#PBS -j eo
#PBS -l nodes=9:ppn=8,walltime=04:00:00
#PBS -q regular
#PBS -A m1871

# Took an hour and 8 minutes

codeLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/"
outputLoc="/Users/Adrian/Research/GlobalSignalInterferometer/gs_mpi/testingFiles2/"
numProcs=4
maxl=2
beam_sig=1.57 
del_bl=0.844 
sqGridSideLen=2
lowerFreq=120 # in MHz
upperFreq=122 # in MHz
deltaFreq=1 # in MHz
variableBeam=0
# 0 is freq-independent primary beam
# 1 is beam size proportional to wavelength, with beam_sig defined to be the size at 150 MHz


# For running on my laptop or Carver
echo "Starting run now..."
date
mpirun -np $numProcs python "$codeLoc/mpi_Q_matrix_grid_mult_fq_3.py" $outputLoc $maxl $beam_sig $del_bl $sqGridSideLen $lowerFreq $upperFreq $deltaFreq $variableBeam
date