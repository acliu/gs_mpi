#! /usr/bin/env python
"""This is a really silly script that's used to make
the necessary directories for the Monte Carlo results""" 

import sys,os,numpy as np

mc_loc = sys.argv[1]
del_bl = float(sys.argv[2])
sqGridSideLen = int(sys.argv[3])
beam_sig = float(sys.argv[4])
lowerFreq = float(sys.argv[5])
upperFreq = float(sys.argv[6])
freqSpace = float(sys.argv[7])
variableBeam = int(sys.argv[8])

fqs = np.arange(lowerFreq,upperFreq+freqSpace,freqSpace)
fqs /= 1000. # Convert from MHz to GHz
numFreqs = fqs.shape[0]

if variableBeam == 0:
    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
elif variableBeam == 1:
    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_lambdaBeam_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)


for fq in fqs:
    if '{0}_fq_{1:.3f}'.format(savekey,fq) not in os.listdir(mc_loc):
        os.mkdir('{0}/{1}_fq_{2:.3f}'.format(mc_loc,savekey,fq))

