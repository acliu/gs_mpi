from mpi4py import MPI 
import sys
import aipy as a, numpy as n
import useful_functions as uf
from scipy import special
import basic_amp_aa_grid_gauss as agg

def get_Q_element_mult_fqs(tx,ty,tz,dOmega,amp,baseline,l,m):
    """
    This returns a vector of Q elements for multiple frequencies, so the 
    vector is the same length as fqs. 
    """
    bx,by,bz = baseline
    # compute spherical harmonic
    Y = n.array(special.sph_harm(m,l,theta,phi)) #using math convention of theta=[0,2pi], phi=[0,pi]    
    #fringe pattern
    tb_grid,fqs_grid = n.meshgrid((bx*tx+by*ty+bz*tz),fqs)
    phs = n.exp(-2j*n.pi*fqs_grid*tb_grid)
    #phs = n.exp(-2j*n.pi*fqs*(bx*tx+by*ty+bz*tz)) # bl in ns, fq in GHz => bl*fq = 1
#    tx.shape = im.uv.shape
    valid = n.logical_not(tx.mask)
#    Y.shape = phs.shape = amp.shape = im.uv.shape
    amp = n.where(valid, amp, n.zeros_like(amp))
    phs = n.where(valid, phs, n.zeros_like(phs))
    Y = n.where(valid, Y, n.zeros_like(Y)) 
    Q_elements = uf.vdot(phs,amp*Y*dOmega) #n.sum(amp*Y*phs*dOmega)
    return Q_elements

def get_dOmega(tx,ty):
    dx = n.zeros_like(tx)
    for ii in range(tx.shape[1]-1):
        dx[:,ii] = n.abs(tx[:,ii]-tx[:,ii+1])
    dx[:,-1] = dx[:,-2]
    dy = n.zeros_like(ty)
    for ii in range(ty.shape[0]-1):
        dy[ii,:] = n.abs(ty[ii,:]-ty[ii+1,:])
    dy[-1,:] = dy[-2,:]
    dOmega = dx*dy/n.sqrt(1-tx*tx-ty*ty)
    return dOmega

#print "everything has imported"

# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1

# define scripts directory location
#save_loc = '/global/scratch2/sd/mpresley/gs_data'
save_loc = '/global/scratch2/sd/acliu/GlobalSignalInterferometer/gs_data'
#save_loc = '/Users/mpresley/Research/Research_Adrian_Aaron/gs_data'

# define parameters related to calculation 
maxl = 10
_,beam_sig,del_bl,num_bl = sys.argv
beam_sig=float(beam_sig); del_bl=float(del_bl);num_bl=int(num_bl)
fqs = n.arange(50,91,2)*0.001

savekey = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}'.format(del_bl,num_bl,beam_sig)

#global tx,ty,tz,dOmega,theta,phi,amp 
im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
dOmega = get_dOmega(tx,ty)
valid = n.logical_not(tx.mask)
tx,ty,tz,dOmega = tx.flatten(),ty.flatten(),tz.flatten(),dOmega.flatten()
theta = n.arctan2(ty,tx) # using math convention of theta=[0,2pi], phi=[0,pi]
phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))
amp = uf.gaussian(beam_sig,n.zeros_like(theta),phi)

baselines = agg.make_pos_array(del_bl,num_bl)

num0,num1 = len(baselines),(maxl+1)*(maxl+1)
print "num baselines = {0}\nnum lms = {1}".format(num0,num1)
lms = n.zeros([num1,2])
ii=0
for ll in range(maxl+1):
	for mm in range(-ll,ll+1):
		lms[ii] = n.array([ll,mm])
		ii+=1
matrix = n.zeros([num0,num1,len(fqs)],dtype=n.complex)
assignment_matrix = n.arange(num0*num1).reshape((num0,num1))

# define parameters related to task-mastering
numToDo = num0*num1
#print "numToDo = ",numToDo
num_sent = 0 # this functions both as a record of how many assignments have 
             # been sent and as a tag marking which matrix entry was calculated
    

# Big running loop
# If I am the master process
if rank==master:
    #print "I am the master! Muahaha!"
    # send out first round of assignments
    for kk in range(num_slaves):
        selectedi, selectedj = n.where(assignment_matrix==num_sent)
        selectedi = selectedi[0]; selectedj = selectedj[0]
        #print "num_sent = ",num_sent
        comm.send(selectedi,dest=kk+1)
        comm.send(selectedj,dest=kk+1)
        #print "i,j = ",selectedi,selectedj," was sent to slave ",kk+1
        num_sent +=1
    #print "Master sent out first round of assignments"
    # listen for results and send out new assignments
    for kk in range(numToDo):
        source,entry = comm.recv(source=MPI.ANY_SOURCE)
        selectedi = comm.recv(source=source)
        selectedj = comm.recv(source=source)
        # stick entry into matrix 
        matrix[selectedi,selectedj,:] = entry
        #print 'Master just received element (i,j) = ',selectedi,selectedj,' from slave ',source
        print 'Have completed {0} of {1}'.format(kk,numToDo)
        # if there are more things to do, send out another assignment
        if num_sent<numToDo:
            selectedi, selectedj = n.where(assignment_matrix==num_sent)
            selectedi = selectedi[0]; selectedj = selectedj[0]
            comm.send(selectedi,dest=source)
            comm.send(selectedj,dest=source)
            #print "Master sent out i,j = ",selectedi, selectedj,' to slave ',source
            num_sent +=1
        else:
            # send a -1 to tell slave that task is complete
            comm.send(-1,dest=source)
            comm.send(-1,dest=source)
            #print "Master sent out the finished i,j to slave ",source
# If I am a slave and there are not more slaves than jobs
elif rank<=numToDo:
    #print "I am slave ",rank
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
            element = get_Q_element_mult_fqs(tx,ty,tz,dOmega,amp,baselines[selectedi],lms[selectedj,0],lms[selectedj,1])
            # send answer back
            comm.send((rank,element),dest=master)
            comm.send(selectedi,dest=master)
            comm.send(selectedj,dest=master)
#            print "Slave ",rank," sent back i,j = ",selectedi,selectedj
comm.Barrier()

if rank==master:
    for ii,fq in enumerate(fqs):
        n.savez_compressed('{0}/Q_matrices/Q_{1}_fq_{2:.3f}'.format(save_loc,savekey,fq),Q=matrix[:,:,ii],baselines=baselines,lms=lms)
    print "The master has saved the matrix."

MPI.Finalize()


