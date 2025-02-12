#print "I exist!"
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
    rb_grid,fqs_grid = n.meshgrid((bix*rx+biy*ry+biz*rz),fqs)
    Gi = amp*n.exp(-2j*n.pi*fqs_grid*rb_grid)*dOmega
    fqs_grid, rb_grid = n.meshgrid((bjx*rx+bjy*ry+bjz*rz),fqs)
    Gj_star = n.conj(amp*n.exp(-2j*n.pi*fqs_grid*rb_grid))*dOmega
    elements = uf.vdot(Gi*Gj_star,Rsq)
    return elements

# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1
#print "defined mpi paramters"

# define file locations
#fits_file_loc = '/global/homes/m/mpresley/scripts/general_files/fits_files/hi1001_32.fits'
fits_file_loc = '/Users/mpresley/soft/gsm/data_50MHz_100MHz/gsm1001_32.fits'
#save_loc = '/global/scratch2/sd/mpresley/gs_data'
save_loc = '/Users/mpresley/Research/Research_Adrian_Aaron/gs_data'

# define parameters related to calculation
fqs = n.arange(50,91,2)*0.001
healmap = a.map.Map(fromfits=fits_file_loc)
global px_array; px_array = n.arange(healmap.npix()) # gets an array of healpix pixel indices
global crd_array; crd_array = n.array(healmap.px2crd(px_array,ncrd=3)) # finds the topocentric coords for each healpix pixel
global Rdata; Rdata = healmap.map.map
global Rsq; Rsq = Rdata*Rdata
global dOmega; dOmega = 4*n.pi/px_array.shape[0]
phi,theta = n.array(healmap.px2crd(px_array,ncrd=2))
#print 'theta max = ',max(theta)
#print 'phi max = ',max(phi)

maxl = 10
_,beam_sig,del_bl,num_bl = sys.argv
beam_sig=float(beam_sig); del_bl=float(del_bl); num_bl=int(num_bl)

savekey = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}'.format(del_bl,num_bl,beam_sig)

amp = uf.gaussian(beam_sig,n.zeros_like(theta),phi)
baselines = agg.make_pos_array(del_bl,num_bl)
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
        print 'Have completed {0} of {1}'.format(kk,numToDo)
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
            print "slave ",rank," acknoledges job completion"
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

