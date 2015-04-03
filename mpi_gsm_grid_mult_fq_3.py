#! /usr/bin/env python

from mpi4py import MPI
import sys
import numpy as n
import useful_functions as uf
import aipy as a
import basic_amp_aa_grid_gauss as agg

#print "Everything imported!"

def compute_element(bli,blj,amp):
    bix,biy,biz = bli; bjx,bjy,bjz = blj
    rx,ry,rz = crd_array
    Gi = amp*n.exp(-2j*n.pi*fq*(bix*rx+biy*ry+biz*rz))*dOmega
    Gj_star = n.conj(amp*n.exp(-2j*n.pi*fq*(bjx*rx+bjy*ry+bjz*rz)))*dOmega
    element = n.sum(Gi*Gj_star*Rsq)
    return element

def compute_element_mult_fqs(bli,blj,amp):
    bix,biy,biz = bli; bjx,bjy,bjz = blj
    rx,ry,rz = crd_array
    rb_grid, fqs_grid = n.meshgrid((bix*rx+biy*ry+biz*rz),fqs)
    Gi = amp*n.exp(-2j*n.pi*fqs_grid*rb_grid)*dOmega
    rb_grid, fqs_grid = n.meshgrid((bjx*rx+bjy*ry+bjz*rz),fqs)
    Gj_star = n.conj(amp*n.exp(-2j*n.pi*fqs_grid*rb_grid))*dOmega
    elements = uf.vdot(Gi*Gj_star,Rsq)
    print "elements shape",elements.shape
    return elements

# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1
#print "defined mpi paramters"

# define file locations
fits_file_loc = sys.argv[1]
#'/Users/mpresley/soft/gsm/data_50MHz_100MHz/gsm1001_32.fits'
save_loc = sys.argv[2]
#'/Users/mpresley/Research/Research_Adrian_Aaron/gs_data'
lowerFreq = float(sys.argv[3])
upperFreq = float(sys.argv[4])
freqSpace = float(sys.argv[5])
maxl = int(sys.argv[6])
beam_sig = float(sys.argv[7])
del_bl = float(sys.argv[8])
sqGridSideLen = int(sys.argv[9])
variableBeam = int(sys.argv[10])

# define parameters related to calculation
fqs = n.arange(lowerFreq,upperFreq+freqSpace,freqSpace)
fqs /= 1000. # Convert from MHz to GHz
if variableBeam == 0:
    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
    beam_sig_fqs = beam_sig * n.ones_like(fqs)
elif variableBeam == 1:
    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_lambdaBeam_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
    #beam_sig_fqs = beam_sig * 0.15 / fqs
    beam_sig_fqs = beam_sig * fqs[0] / fqs

healmap = a.map.Map(fromfits=fits_file_loc)
global px_array; px_array = n.arange(healmap.npix()) # gets an array of healpix pixel indices
global crd_array; crd_array = n.array(healmap.px2crd(px_array,ncrd=3)) # finds the topocentric coords for each healpix pixel
global Rdata; Rdata = healmap.map.map
global Rsq; Rsq = Rdata*Rdata
global dOmega; dOmega = 4*n.pi/px_array.shape[0]
phi,theta = n.array(healmap.px2crd(px_array,ncrd=2))
#print 'theta max = ',max(theta)
#print 'phi max = ',max(phi)

amp = n.zeros((fqs.shape[0],len(phi)))
for i,beamSize in enumerate(beam_sig_fqs):
    #amp[i,:] = uf.gaussian(beamSize,n.zeros_like(theta),phi)
    amp[i,:] = uf.cosine_gaussian(beamSize,phi)
#baselines = agg.make_pos_array(del_bl,sqGridSideLen)
baselines = agg.make_uhp_bls(del_bl,sqGridSideLen)

#print "defined calculation parameters"

# define matrix to be calculated
num = len(baselines)
matrix = n.zeros([num,num,len(fqs)],dtype=n.complex)
# define parameters related to task-mastering
numToDo = num*(num+1)/2
#print 'numToDo = ',numToDo
assn_inds = []
for ii in range(num+1):
    for jj in range(ii+1):
        assn_inds.append((ii,jj))
num_sent = 0 # this functions both as a record of how many assignments have 
            # been sent and as a tag marking which matrix entry was calculated
#print "just before the big if statement"

# Big running loop
# If I am the master process
if rank==master:
#    print "I am the master! Muahaha!"
    # send out first round of assignments
    for kk in range(num_slaves):
        selectedi, selectedj = assn_inds[kk]
#        print "num_sent = ",num_sent
        comm.send(selectedi,dest=kk+1)
        comm.send(selectedj,dest=kk+1)
#        print "i,j = ",selectedi,selectedj," was sent to slave ",kk+1
        num_sent +=1
#    print "Master sent out first round of assignments"
    # listen for results and send out new assignments
    for kk in range(numToDo):
        source,entry = comm.recv(source=MPI.ANY_SOURCE)
        selectedi = comm.recv(source=source)
        selectedj = comm.recv(source=source)
        # stick entry into matrix 
        matrix[selectedi,selectedj,:] = entry
        matrix[selectedj,selectedi,:] = n.conj(entry)
        print 'Master just received element (i,j) = ',selectedi,selectedj,' from slave ',source
        print 'Have completed {0} of {1}'.format(kk+1,numToDo)
        # if there are more things to do, send out another assignment
        if num_sent<numToDo:
            selectedi, selectedj = assn_inds[num_sent]
            comm.send(selectedi,dest=source)
            comm.send(selectedj,dest=source)
#            print "Master sent out i,j = ",selectedi, selectedj,' to slave ',source
            num_sent +=1
        else:
            # send a -1 to tell slave that task is complete
            comm.send(-1,dest=source)
            comm.send(-1,dest=source)
#            print "Master sent out the finished i,j to slave ",source
# If I am a slave and there are not more slaves than jobs
elif rank<=numToDo:
#    print "I am slave ",rank
    complete = False
    while not complete:
        # Get assignment
        selectedi = comm.recv(source=master)
        selectedj = comm.recv(source=master)
#        print "slave ",rank," just recieved i,j = ",selectedi,selectedj
        if selectedi==-1:
            # if there are no more jobs
            complete=True
            print "slave ",rank," acknowledges job completion"
        else:
            # compute the matrix element
            bli = baselines[selectedi,:]
            blj = baselines[selectedj,:]
            element = compute_element_mult_fqs(bli,blj,amp)
            # send answer back
            comm.send((rank,element),dest=master)
            comm.send(selectedi,dest=master)
            comm.send(selectedj,dest=master)
#            print "Slave ",rank," sent back i,j = ",selectedi,selectedj
comm.Barrier()

if rank==master:
    for ii,fq in enumerate(fqs):
        n.savez_compressed('{0}/gsm_matrices/gsm_{1}_fq_{2:.3f}'.format(save_loc,savekey,fq),matrix=matrix[:,:,ii],baselines=baselines)
    print "The master has saved the matrix."

MPI.Finalize()

