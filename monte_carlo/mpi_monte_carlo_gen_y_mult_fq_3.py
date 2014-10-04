#! /usr/bin/env python

from mpi4py import MPI 
import sys, os
import numpy as n, aipy as a
import useful_functions as uf
import basic_amp_aa_grid_gauss as agg

def haslam_extrap(chunkLen,randos):
    """
    Returns an array with realizations of the real sky.
    The array measures numRealizations x numFreqs x numPixels
    """
    # Scale the random numbers given the variance of each perturbation
    for i in range(chunkLen):
        for j in range(numPerts):
            randos[i,j] *= pertVar[j]
    #if rank==3:
    #    print "new randoms",randos[-1]
    #    print "freqRatioPow",freqRatioPow
    #    print "haslamBasicExtrap",haslamBasicExtrap
    simulatedSkies = n.exp(n.einsum('ijk,mj',randos,freqRatioPow))
    simulatedSkies = n.transpose(simulatedSkies,(0,2,1))
    simulatedSkies *= haslamBasicExtrap
    #if rank==3:
    #    print "simulation",simulatedSkies[-1]
    return simulatedSkies # MCs x freqs x pixels

def extractSampleSkies(chunkNum,simSkies):
    """
    Outputs a sample sky and a summed sky (over spatial directions).
    """
    sampleSky = simSkies[0]
    outMap = a.map.Map()
    for i,sample, in enumerate(sampleSky):
        outMap.set_map(sample)
        outMap.to_fits('{0}/{1}_fq_{2:.3f}/{1}_fq_{2:.3f}_chunk_{3}_sample.fits'\
                         .format(mc_loc,savekey,fqs[i],chunkNum),clobber=True)
    
    monopoleSky = n.sum(simSkies,axis=2) / float(npix) # MCs x freqs
    n.save('{0}/spatialMean_{1}_chunk_{2}.npy'.format(mc_loc,savekey,chunkNum),monopoleSky)
    return None

def simulateSingleDipole(chunkNum,simSkies):
    localNumMCs = simSkies.shape[0]
    numFreqs = simSkies.shape[1]
    singleDipoleResult = n.zeros((localNumMCs,numFreqs))
    for i,sky in enumerate(simSkies):
        for j,(beam,sk) in enumerate(zip(primaryBeam,sky)):
            singleDipoleResult[i,j] = n.dot(beam,sk) # MCs x freqs
    for i in range(localNumMCs):
        singleDipoleResult[i] /= integratedBeams

    n.save('{0}/singleDipoleResult_{1}_chunk_{2}.npy'.format(mc_loc,savekey,chunkNum),singleDipoleResult)
    return None

def simulate_measurement(chunkLen,simSkies):
    yVects = n.zeros((chunkLen,numFreqs,numBl),dtype=complex)
    for i in range(chunkLen):
        for j in range(numFreqs):
            yVects[i,j] = n.einsum('ab,b',Gmatrices[j],simSkies[i,j])
    return yVects # MCs x freqs x bls

# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1

mc_loc = sys.argv[1]
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
numMCs = int(sys.argv[14])
saveInterval = int(sys.argv[15])

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
pertVar = n.delete(pertVar,0,0) # Remove the unperturbed spectral index because we don't need it
pertVar = n.sqrt(pertVar) # Turn the variance into a standard deviation

# The idea behind this code is that the full list of Monte Carlos is split
# into chunks given by the saveInterval.  Each child process gets a chunk
# to do at a time, and then just writes directly to disk rather than sending
# everything back to the master.  It's best to have numMC / saveInterval >>
# nprocs

numChunks = numMCs / saveInterval
if numMCs%saveInterval != 0:
    numChunks += 1
assignmentArray = n.zeros((numChunks,2),dtype=int)
for i in range(numChunks-1):
    assignmentArray[i,0] = i * saveInterval
    assignmentArray[i,1] = (i+1) * saveInterval - 1
assignmentArray[-1,0] = (numChunks-1) * saveInterval
assignmentArray[-1,1] = numMCs - 1

# If I am the master process
if rank==master:
    print "I am the master! Muahaha!"
    num_sent = 0
    # send out first round of assignments
    for kk in range(min(numChunks,num_slaves)):
        selectedChunkNum = num_sent
        print "num_sent = ",num_sent
        comm.send(selectedChunkNum,dest=kk+1)
        #print "i = ",selectedi," was sent to slave ",kk+1
        num_sent +=1
    print "Master sent out first round of assignments"
    # listen for results and send out new assignments
    for kk in range(numChunks):
        source,selectedChunkNum = comm.recv(source=MPI.ANY_SOURCE)
        print "Master just got chunkNum=",selectedChunkNum
        # if there are more things to do, send out another assignment
        if num_sent<numChunks:
            selectedChunkNum = num_sent
            comm.send(selectedChunkNum,dest=source)
            print "Master sent out i = ",selectedChunkNum,' to slave ',source
            num_sent +=1
        else:
            # send a -1 to tell slave that task is complete
            comm.send(-1,dest=source)
            print "Master sent out the finished i to slave ",source
elif rank<=numChunks:
    print "I am slave ",rank
    complete = False
    while not complete:
        # Get assignment
        selectedChunkNum = comm.recv(source=master)
        #print "slave ",rank," just recieved i = ",selectedi
        if selectedChunkNum == -1:
            # if there are no more jobs
            complete=True
            print "slave ",rank," acknowledges job completion"
        else:
            # compute the matrix element
            print "slave working on chunkNum",selectedChunkNum
            # send answer back
            startIndex,endIndex = assignmentArray[selectedChunkNum]
            chunkLen = endIndex - startIndex + 1
            n.random.seed(startIndex)
            randNums = n.random.normal(size=chunkLen*numPerts*npix)
            randNums = randNums.reshape((chunkLen,numPerts,npix))
            simSkies = haslam_extrap(chunkLen,randNums)
            extractSampleSkies(selectedChunkNum,simSkies)
            simulateSingleDipole(selectedChunkNum,simSkies)
            simulated_yVects = simulate_measurement(chunkLen,simSkies)
            for i,freq in enumerate(fqs):
                n.savez_compressed('{0}/{1}_fq_{2:.3f}/mc_{1}_chunk_{3}'.\
                                    format(mc_loc,savekey,freq,selectedChunkNum),\
                                    matrix=simulated_yVects[:,i,:])
            comm.send((rank,selectedChunkNum),dest=master)
            print "Slave ",rank," sent back i = ",selectedChunkNum
comm.Barrier()
MPI.Finalize()

