#! /usr/bin/env python
"""This script takes in the Monte Carlo 'chunks' and consolidates
them into single files.""" 

import sys,os,numpy as np

mc_loc = sys.argv[1]
del_bl = float(sys.argv[2])
sqGridSideLen = int(sys.argv[3])
beam_sig = float(sys.argv[4])
lowerFreq = float(sys.argv[5])
upperFreq = float(sys.argv[6])
freqSpace = float(sys.argv[7])
variableBeam = int(sys.argv[8])
numMCs = int(sys.argv[9])
saveInterval = int(sys.argv[10])

fqs = np.arange(lowerFreq,upperFreq+freqSpace,freqSpace)
fqs /= 1000. # Convert from MHz to GHz
numFreqs = fqs.shape[0]

numChunks = numMCs / saveInterval
if numMCs%saveInterval != 0:
    numChunks += 1

if variableBeam == 0:
    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
elif variableBeam == 1:
    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_lambdaBeam_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)

for i,freq in enumerate(fqs):
    fullMCs = None
    for selectedChunkNum in range(numChunks):
        temp = np.load('{0}/{1}_fq_{2:.3f}/mc_{1}_chunk_{3}.npz'.format(mc_loc,savekey,freq,selectedChunkNum))
        if fullMCs == None:
            fullMCs = temp['matrix']
        else:
            fullMCs = np.vstack((fullMCs,temp['matrix']))
        temp.close()
    np.savez_compressed('{0}/{1}_fq_{2:.3f}/mc_{1}_fq_{2:.3f}_allMCs'.format(mc_loc,savekey,freq),matrix=fullMCs)

monopoleMCs = None
for selectedChunkNum in range(numChunks):
    temp = np.load('{0}/spatialMean_{1}_chunk_{2}.npy'.format(mc_loc,savekey,selectedChunkNum))
    if monopoleMCs == None:
        monopoleMCs = temp
    else:
        monopoleMCs = np.vstack((monopoleMCs,temp))

ensembleAvMonopole = np.mean(monopoleMCs,axis=0)
np.save('{0}/spatialMean_{1}_ensembleAv.npy'.format(mc_loc,savekey),ensembleAvMonopole)

covariance = np.zeros((numFreqs,numFreqs))
for realization in monopoleMCs:
    covariance += np.outer(realization,realization)
covariance /= float(numMCs)
covariance -= np.outer(ensembleAvMonopole,ensembleAvMonopole)

np.save('{0}/covariance_{1}.npy'.format(mc_loc,savekey),covariance)