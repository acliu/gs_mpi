#! /usr/bin/env python
"""Does a single realization of the monte carlo, extrapolating
from Haslam, outputting the maps for later use as well as
generating the simulated visibilities."""

import sys
import numpy as n, aipy as a
import useful_functions as uf
import basic_amp_aa_grid_gauss as agg

def haslam_extrap(randos):
    """
    Returns an array with realizations of the real sky.
    The array measures numFreqs x numPixels
    """

    # Scale the random numbers given the variance of each perturbation
    for j in range(numPerts):
        randos[j] *= pertVar[j]
    
    print "here are the randos"
    for i in range(randos.shape[0]):
        print randos[i]

    print "here is the freqRatioPow before"
    for i in range(freqRatioPow.shape[0]):
        print freqRatioPow[i]

    simulatedSkies = n.exp(n.einsum('jk,mj',randos,freqRatioPow))

    print "here is the freqRatioPow after"
    for i in range(simulatedSkies.shape[0]):
        print simulatedSkies[i]

    simulatedSkies = simulatedSkies.T
    simulatedSkies *= haslamBasicExtrap

    return simulatedSkies # freqs x pixels

def simulateSingleDipole(sky):
    numFreqs = sky.shape[0]
    singleDipoleResult = n.zeros(numFreqs)
    for j,(beam,sk) in enumerate(zip(primaryBeam,sky)):
        singleDipoleResult[j] = n.dot(beam,sk) # freqs
    singleDipoleResult /= integratedBeams

    n.save('{0}/singleRealization_singleDipoleResult_{1}.npy'.format(output_loc,savekey),singleDipoleResult)
    return None

def simulate_measurement(sky):
    yVects = n.zeros((numFreqs,numBl),dtype=complex)
    for j in range(numFreqs):
        yVects[j] = n.einsum('ab,b',Gmatrices[j],sky[j])
    return yVects # freqs x bls

output_loc = sys.argv[1]
hasmap_loc = sys.argv[2]
Gmatrix_loc = sys.argv[3]
nside = int(sys.argv[4])
del_bl = float(sys.argv[5])
sqGridSideLen = int(sys.argv[6])
beam_sig = float(sys.argv[7])
variableBeam = int(sys.argv[8])
lowerFreq = float(sys.argv[9])
upperFreq = float(sys.argv[10])
freqSpace = float(sys.argv[11])
numPerts = int(sys.argv[12])
pertFile = sys.argv[13]

fqs = n.arange(lowerFreq,upperFreq+freqSpace,freqSpace)
fqs /= 1000. # Convert from MHz to GHz
numFreqs = fqs.shape[0]
npix = 12 * nside * nside
baselines = agg.make_uhp_bls(del_bl,sqGridSideLen)
numBl = len(baselines)

# This code does not use del_bl, sqGridSideLen, or beam_sig for
# anything other than filename definitions
if variableBeam == 0:
    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
elif variableBeam == 1:
    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_lambdaBeam_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
Gmatrices = n.zeros((numFreqs,numBl,npix),dtype=complex)
for i,freq in enumerate(fqs):
    tempMatrix = n.load('{0}/G_matrices/G_{1}_fq_{2:.3f}.npz'.format(Gmatrix_loc,savekey,freq))
    Gmatrices[i] = tempMatrix['matrix'].copy() # All processes read in the G matrices
    tempMatrix.close()

# Create primary beams for everyone to use
beamMap = a.map.Map(nside)
directionVects_xyz = beamMap.map.px2crd(n.array([i for i in range(npix)]))
directionVects_xyz = n.array(directionVects_xyz).T
directionVects_thetas = beamMap.map.px2crd(n.array([i for i in range(npix)]),ncrd=2)[0]
    
if variableBeam == 0:
    beam_sig_fqs = beam_sig * n.ones_like(fqs)
elif variableBeam == 1:
    beam_sig_fqs = beam_sig * 0.15 / fqs
   
primaryBeam = n.zeros((fqs.shape[0],npix))
for i,beamSize in enumerate(beam_sig_fqs):
    #primaryBeam[i,:] = uf.gaussian(beamSize,n.zeros(npix),directionVects_thetas)
    primaryBeam[i,:] = uf.cosine_gaussian(beamSize,directionVects_thetas)

integratedBeams = n.sum(primaryBeam,axis=1)

# Everyone reads in the Haslam map and forms the extrapolation to lower frequencies
pertVar = n.loadtxt(pertFile) # First component of this is the unperturbed spectral index
nu0=0.408 # 408 MHz
hasmap = a.map.Map(fromfits=hasmap_loc)
hasmapASCII = hasmap.map.get_map()
freqRatio = fqs / nu0
freqRatioLog = n.log(freqRatio.copy())
freqRatioPow = freqRatioLog.copy()
if numPerts > 1:
    for i in range(2,numPerts+1):
        freqRatioPow = n.vstack((freqRatioPow,freqRatioLog**i))
    freqRatioPow = freqRatioPow.T
haslamBasicExtrap = n.outer(freqRatio**pertVar[0],hasmapASCII)
n.save('{0}/seriously_{1}.npy'.format(output_loc,savekey),freqRatio**pertVar[0])
n.save('{0}/nullSkyTest_{1}.npy'.format(output_loc,savekey),haslamBasicExtrap)
pertVar = n.delete(pertVar,0,0) # Remove the unperturbed spectral index because we don't need it
pertVar = n.sqrt(pertVar) # Turn the variance into a standard deviation

n.random.seed(5678)
randNums = n.random.normal(size=numPerts*npix)
randNums = randNums.reshape((numPerts,npix))
simSkies = haslam_extrap(randNums)
n.save('{0}/singleRealization_skyMaps_{1}.npy'.format(output_loc,savekey),simSkies)
simulateSingleDipole(simSkies)
simulated_yVects = simulate_measurement(simSkies)
for i,freq in enumerate(fqs):
    n.savez_compressed('{0}/{1}_fq_{2:.3f}/singleRealization_y_{1}'.\
                    format(output_loc,savekey,freq),\
                    matrix=simulated_yVects[i,:])
