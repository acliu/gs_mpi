#! /usr/bin/env python

import numpy as np, sys
import basic_amp_aa_grid_gauss as agg

def vect2sq(n,vect):
    sq = np.zeros((n,n),dtype=complex)
    for i in range(n):
        for j in range(n):
            if i >= j:
                sq[i,j] = vect[j+i*(i+1)/2]
            else:
                sq[i,j] = np.conj(vect[i+j*(j+1)/2])
    return sq

if __name__=='__main__':
    nside = int(sys.argv[1])
    npix = 12 * nside * nside
    KvectFname = sys.argv[2]
    GmatrixLoc = sys.argv[3]
    del_bl = float(sys.argv[4])
    sqGridSideLen = int(sys.argv[5])
    beam_sig = float(sys.argv[6])
    variableBeam = int(sys.argv[7])
    lowerFreq = float(sys.argv[8])
    upperFreq = float(sys.argv[9])
    freqSpace = float(sys.argv[10])
    save_loc = sys.argv[11]

    fqs = np.arange(lowerFreq,upperFreq+freqSpace,freqSpace)
    fqs /= 1000. # Convert from MHz to GHz
    numFreqs = fqs.shape[0]
    baselines = agg.make_uhp_bls(del_bl,sqGridSideLen)
    numBl = len(baselines)

    if variableBeam == 0:
        savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
    elif variableBeam == 1:
        savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_lambdaBeam_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
    Gmatrices = np.zeros((numFreqs,numBl,npix),dtype=complex)
    for i,freq in enumerate(fqs):
        tempMatrix = np.load('{0}/G_matrices/G_{1}_fq_{2:.3f}.npz'.format(GmatrixLoc,savekey,freq))
        Gmatrices[i] = tempMatrix['matrix'] # All processes read in the G matrices
    Kvect = np.load(KvectFname)
    Kmatrix = vect2sq(npix,Kvect)

    for ii,fq in enumerate(fqs):
        Gmatrix = Gmatrices[ii]
        GKmatrix = np.einsum('ij,jk',Gmatrix,Kmatrix)
        GKGdagger = np.einsum('ij,jk',GKmatrix,np.conj(Gmatrix.T))
        np.savez_compressed('{0}/gsm_matrices/improvedNfg_{1}_fq_{2:.3f}'.format(save_loc,savekey,fq),matrix=GKGdagger,baselines=baselines)
