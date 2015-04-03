#! /usr/bin/env python

import aipy as ap, numpy as np, sys
import useful_functions as uf
import basic_amp_aa_grid_gauss as agg


if __name__=='__main__': 
    nside = int(sys.argv[1])
    npix = 12 * nside * nside
    inputMapFname = sys.argv[2]
    outputMapFname = sys.argv[3]
    beam_sig = 1.57 * 0.15 / 0.12 * 0.4

    beamMap = ap.map.Map(nside)
    directionVects_xyz = beamMap.map.px2crd(np.array([i for i in range(npix)]))
    directionVects_xyz = np.array(directionVects_xyz).T
    directionVects_thetas = beamMap.map.px2crd(np.array([i for i in range(npix)]),ncrd=2)[0]

    primaryBeam = uf.gaussian(beam_sig,np.zeros(npix),directionVects_thetas)
    totalBeamArea = np.sum(primaryBeam)

    inputMap = ap.map.Map(fromfits=inputMapFname)
    inputASCII = inputMap.map.get_map()
    simpleAverage = np.sum(inputASCII) / npix
    inputASCII *= primaryBeam
    inputMap.map.set_map(inputASCII)
    inputMap.to_fits(outputMapFname)

    integratedMonopoleResponse = np.sum(inputASCII)
    singleDipoleExtraction = integratedMonopoleResponse / totalBeamArea

    print "simple extraction gives:",singleDipoleExtraction
    print "simple averaging gives:",simpleAverage

    # lowerFreq = float(sys.argv[2])
    # upperFreq = float(sys.argv[3])
    # freqSpace = float(sys.argv[4])
    # beam_sig = float(sys.argv[5])
    # del_bl = float(sys.argv[6])
    # sqGridSideLen = int(sys.argv[7])
    # variableBeam = int(sys.argv[8])
    # save_loc = sys.argv[9]

    # if variableBeam == 0:
    #     savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
    # elif variableBeam == 1:
    #     savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_lambdaBeam_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)

    # #baselines = agg.make_pos_array(del_bl,sqGridSideLen)
    # baselines = agg.make_uhp_bls(del_bl,sqGridSideLen)
    # numBl = len(baselines)
    # dOmega = 4 * np.pi / npix

    # beamMap = ap.map.Map(nside)
    # directionVects_xyz = beamMap.map.px2crd(np.array([i for i in range(npix)]))
    # directionVects_xyz = np.array(directionVects_xyz).T
    # directionVects_thetas = beamMap.map.px2crd(np.array([i for i in range(npix)]),ncrd=2)[0]
    
    # fqs = np.arange(lowerFreq,upperFreq+freqSpace,freqSpace)
    # fqs /= 1000. # Convert from MHz to GHz

    # if variableBeam == 0:
    #     beam_sig_fqs = beam_sig * np.ones_like(fqs)
    # elif variableBeam == 1:
    #     beam_sig_fqs = beam_sig * 0.15 / fqs
    
    # primaryBeam = np.zeros((fqs.shape[0],npix))
    # for i,beamSize in enumerate(beam_sig_fqs):
    #     primaryBeam[i,:] = uf.gaussian(beamSize,np.zeros(npix),directionVects_thetas)

    # for k,freq in enumerate(fqs):
    #     print "Doing frequency number",k,"which is at",freq,"GHz"
    #     Gmatrix = np.zeros((numBl,npix),dtype=complex)
    #     for j,(nhat,beamVal) in enumerate(zip(directionVects_xyz,primaryBeam[k])):
    #         for i,bl in enumerate(baselines):
    #             Gmatrix[i,j] = np.exp(-2j*np.pi*freq*np.dot(nhat,bl))
    #         Gmatrix[:,j] *= beamVal
    #     Gmatrix *= dOmega
    #     np.savez_compressed('{0}/G_matrices/G_{1}_fq_{2:.3f}'.format(save_loc,savekey,freq),matrix=Gmatrix)