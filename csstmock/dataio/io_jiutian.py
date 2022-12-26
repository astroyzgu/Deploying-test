'''
This "io_jiutian" module is writen for extracting data from Jiutian simulation. 

collect4csstmock: perpare data for the construction of mock catalog.  

read_groups:  extract halo/subhalo informations of one snapshot from all subfind group files.

read_hbt: read subhalo information of one snapshot from all SubSnap files. 

nexthaloid: determine the host halo id in next snapshot.

Wirtten by Qingyang Li, Yizhou Gu (2022/11) 
''' 

import glob
import multiprocessing as mp
import numpy as np
import time
import os

# __ALL__ = ['collect4csstmock', 'read_groups', 'read_hbt', 'nexthaloid']

def read_hbt_head(filename):
    '''
    Read some headers from hbt files. 

    Note: 

    Cosmology in all files are same at same snapshot.
    '''
    HEAD = {}
    data = h5py.File(filename,'r')
    HEAD['HubbleParam'] = data['Cosmology/HubbleParam'][0]
    HEAD['OmegaLambda0'] = data['Cosmology/OmegaLambda0'][0]
    HEAD['OmegaM0'] = data['Cosmology/OmegaM0'][0]
    HEAD['ScaleFactor'] = data['Cosmology/ScaleFactor'][0]
    HEAD['NumberOfSubhalosInAllFiles'] = data['NumberOfSubhalosInAllFiles'][0]
    HEAD['NumberOfFiles'] = data['NumberOfFiles'][0]

    #we supplyment information here
    HEAD['NumberOfSubhalosThisFile'] = data['Subhalos'].fields('TrackId')[:].shape[0]
    data.close()

    return HEAD


def read_hbt_subhalos(filename, blocks):
    '''
    read subhalo information from one SubSnap file

    Wirtten by Qingyang Li

    Parameters
    ----------
    filename: path of file
    blocks: name of blocks

    Returns
    -------
    Nsub: subhalo number
    DATA: dict, datasets of blocks
    DATATYPE: dict, data type of blocks 
    '''
    
    import numpy as np
    import h5py

    #all names of blocks
    names = ['TrackId', 'Nbound', 'Mbound', 'HostHaloId', 'Rank', 'Depth', 'LastMaxMass', 'SnapshotIndexOfLastMaxMass', 'SnapshotIndexOfLastIsolation', 'SnapshotIndexOfBirth', 'SnapshotIndexOfDeath', 'SnapshotIndexOfSink', 'RmaxComoving', 'VmaxPhysical', 'LastMaxVmaxPhysical', 'SnapshotIndexOfLastMaxVmax', 'R2SigmaComoving', 'RHalfComoving', 'BoundR200CritComoving', 'BoundM200Crit', 'SpecificSelfPotentialEnergy', 'SpecificSelfKineticEnergy', 'SpecificAngularMomentum', 'InertialEigenVector', 'InertialEigenVectorWeighted', 'InertialTensor', 'InertialTensorWeighted', 'ComovingAveragePosition', 'PhysicalAverageVelocity', 'ComovingMostBoundPosition', 'PhysicalMostBoundVelocity', 'MostBoundParticleId', 'SinkTrackId']
    types = ['<i8', '<i8', '<f4', '<i8', '<i8', '<i4', '<f4', '<i4', '<i4', '<i4', '<i4', '<i4', '<f4', '<f4', '<f4', '<i4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<f4', '<i8', '<i8']

    data = h5py.File(filename,'r')
    nbound= data['Subhalos'].fields('Nbound')[:]
    idx = np.where(nbound > 1)[0]
    Nsub = idx.shape[0]
    
    DATA = {}
    DATATYPE = {}
    for block in blocks:
        readdata = data['Subhalos'].fields(block)[:]
        readdata = readdata[idx] #select Nbound != 1
        nblock   = names.index(block)
        datatype = types[nblock]
        # datatype = data['Subhalos'].dtype[nblock]
        # datatype = data['Subhalos'].dtype.descr[nblock][1]
        # print('Read subhalo blocks: %s,' %block + ' dtype:', datatype)

        DATA[block] = np.array(readdata, dtype = datatype)
        DATATYPE[block] = datatype
    data.close()
    return Nsub, DATA, DATATYPE

def read_hbt_subhalos_split(filenames, blocks):
    '''
    The worker of read_hbt_subhalos for parallel 
    '''
    DATA = {}; DATATYPE = {} 
    for block in blocks:  DATA[block] = []
    
    Nsub = []; 
    for filename in filenames:
        #print(filename)
        Nsub_, DATA_, DATATYPE_ = read_hbt_subhalos(filename, blocks) 
        for block in blocks: DATA[block].append(DATA_[block])
        Nsub.extend([Nsub_])
    DATATYPE = DATATYPE_
    Nsub = np.sum(Nsub, dtype = np.int64)
    blockin3D = ['ComovingAveragePosition', 'PhysicalAverageVelocity', 'ComovingMostBoundPosition', 'PhysicalMostBoundVelocity'] 
    for block in blocks: 
        if block in blockin3D: 
            DATA[block] = np.vstack(DATA[block])
        else: 
            DATA[block] = np.hstack(DATA[block])
    return Nsub, DATA, DATATYPE

def read_hbt(snapnum, blocks, basedir_hbt = '/home/cossim/Jiutian/M1000/hbt/'): 
    '''
    read subhalo information from all SubSnap file with multiprocessing
    Wirtten by Qingyang Li

    Parameters: 
    ----------
    snapnum: int
          the number of snapshot
    blocks:  list
          the name of blocks
    basedir_groups: str 
          file path for subfind data (FoF)

    Returns:
    -------
    DATAALL: dict  
          datasets of blocks. for example, Nbound = DATAALL['Nbound'][:]
    DATATYPE: dict 
          data type of blocks


    Note also
    ---------
    Subhalos with Nbound<=1 are not included. 
    ''' 
    
    #initialize data
    DATAALL = {}
    for block in blocks: DATAALL[block] = []
    filedir   = basedir_hbt + str(snapnum).zfill(3) + \
                '/SubSnap_'+ str(snapnum).zfill(3)  + '.'
    filenames = glob.glob(filedir + '*.hdf5') 
    division  = len(filenames) # the divided  number
    filenames = [filedir + '%s.hdf5'%d for d in range(division) ] 
    
    size        = mp.cpu_count() 
    if division == 0:
        print('No found in ' + filedir)
        exit()
    else:
        print('Reading %s hbt divisions of %sth snapshot with %s threads'%(division, snapnum, size), ) 

    #split task & start parallel
    task_splits = np.array_split( filenames, size)
    pool = mp.Pool(size) 
    res  = []
    for ii in range(size):
        r = pool.apply_async(read_hbt_subhalos_split, args = (task_splits[ii], blocks) ) 
        res.append(r)
    pool.close() 
    pool.join()
    NSUB = []
    for r in res:
        Nsub_, DATA, DATATYPE = r.get() 
        for block in blocks: DATAALL[block].append(DATA[block])
        NSUB.append(Nsub_) 
    NSUB = np.sum(NSUB, dtype = np.int64)
    blockin3D = ['ComovingAveragePosition', 'PhysicalAverageVelocity', 'ComovingMostBoundPosition', 'PhysicalMostBoundVelocity'] 
    for block in blocks:  
        if block in blockin3D: 
            DATAALL[block] = np.vstack(DATAALL[block])
        else: 
            DATAALL[block] = np.hstack(DATAALL[block]) 
    return DATAALL, DATATYPE

class subfind_catalog:
  '''
  code for reading Subfind's subhalo_tab files
  '''
  def __init__(self, curfile): 
 
    #self.filebase = basedir + "/groups_" + str(snapnum).zfill(3) + "/subhalo_tab_" + str(snapnum).zfill(3) + "."
 
    #print()
    #print("reading subfind catalog for snapshot",snapnum,"of",basedir)
 
    have_veldisp = True
 
    #curfile = self.filebase + str(filenum)
    
    if (not os.path.exists(curfile)):
      print("file not found:", curfile)
      sys.exit()
    
    f = open(curfile,'rb')
    
    ngroups = np.fromfile(f, dtype=np.uint32, count=1)[0]    # Number of groups within this file
    totngroups = np.fromfile(f, dtype=np.uint64, count=1)[0] # Total number of groups for this snapshot.
    nids = np.fromfile(f, dtype=np.uint32, count=1)[0]       
    totnids = np.fromfile(f, dtype=np.uint64, count=1)[0]
    ntask   = np.fromfile(f, dtype=np.uint32, count=1)[0]
    nsubs   = np.fromfile(f, dtype=np.uint32, count=1)[0]
    totnsubs = np.fromfile(f, dtype=np.uint64, count=1)[0]
    
    self.ngroups= ngroups # Number of    groups within this file chunk.
    self.nids   = nids  # ???? print(self.nids)  
    self.nfiles = ntask # Total number of file chunks the group catalog is split between.
    self.nsubs  = nsubs # Number of subgroups within this file chunk.

    self.totngrp = totngroups # Total number of    groups for this snapshot.
    self.totnsub = totnsubs   # Total number of subgroups for this snapshot.


    self.group_len = np.empty(ngroups, dtype=np.uint32)
    self.group_offset = np.empty(ngroups, dtype=np.uint32)
    self.group_nr = np.empty(ngroups, dtype=np.uint64)
    self.group_cm = np.empty(ngroups, dtype=np.dtype((np.float32,3)))
    self.group_vel = np.empty(ngroups, dtype=np.dtype((np.float32,3)))
    self.group_pos = np.empty(ngroups, dtype=np.dtype((np.float32,3)))
    self.group_m_mean200 = np.empty(ngroups, dtype=np.float32)
    self.group_m_crit200 = np.empty(ngroups, dtype=np.float32)
    self.group_m_tophat200 = np.empty(ngroups, dtype=np.float32)
    self.group_veldisp = np.empty(ngroups, dtype=np.float32)
    if have_veldisp:
      self.group_veldisp_mean200 = np.empty(ngroups, dtype=np.float32)
      self.group_veldisp_crit200 = np.empty(ngroups, dtype=np.float32)
      self.group_veldisp_tophat200 = np.empty(ngroups, dtype=np.float32)
    self.group_nsubs = np.empty(ngroups, dtype=np.uint32)
    self.group_firstsub = np.empty(ngroups, dtype=np.uint32)
    
    if nsubs > 0: 
      self.sub_len = np.empty(nsubs, dtype=np.uint32)
      self.sub_offset = np.empty(nsubs, dtype=np.uint32)
      self.sub_grnr = np.empty(nsubs, dtype=np.uint64)
      self.sub_nr = np.empty(nsubs, dtype=np.uint64)
      self.sub_pos = np.empty(nsubs, dtype=np.dtype((np.float32,3)))
      self.sub_vel = np.empty(nsubs, dtype=np.dtype((np.float32,3)))
      self.sub_cm = np.empty(nsubs, dtype=np.dtype((np.float32,3)))
      self.sub_spin = np.empty(nsubs, dtype=np.dtype((np.float32,3)))
      self.sub_veldisp = np.empty(nsubs, dtype=np.float32)
      self.sub_vmax = np.empty(nsubs, dtype=np.float32)
      self.sub_vmaxrad = np.empty(nsubs, dtype=np.float32)
      self.sub_halfmassrad = np.empty(nsubs, dtype=np.float32)
      #self.sub_shape = np.empty(nsubs, dtype=np.dtype((np.float32,6)))
      self.sub_ebind = np.empty(nsubs, dtype=np.float32)
      self.sub_pot = np.empty(nsubs, dtype=np.float32)
      #self.sub_profile = np.empty(nsubs, dtype=np.dtype((np.float32,9)))
      self.sub_parent = np.empty(nsubs, dtype=np.uint32)
      self.sub_idbm = np.empty(nsubs, dtype=np.uint64)
    #--------------------------------------------------------------------
    self.group_len = np.fromfile(f, dtype=np.uint32, count=ngroups)
    # group_len - Integer counter of the total number of particles/cells of all types in this group.
    self.group_offset = np.fromfile(f, dtype=np.uint32, count=ngroups)
    self.group_nr = np.fromfile(f, dtype=np.uint64, count=ngroups)

    self.group_cm = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=ngroups)
    self.group_vel = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=ngroups)
    self.group_pos = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=ngroups)
    # group_cm  (N,3) - Center of mass of the group (mass weighted relative coordinates of all particles/cells in the group)
    # group_vel (N,3) - Velocity of the group, The peculiar velocity is obtained by multiplying this value by 1/a. 
    # group_pos (N,3) - Spatial position within the periodic box (of the particle with the minimum gravitational potential energy)
    self.group_m_mean200 = np.fromfile(f, dtype=np.float32, count=ngroups)
    self.group_m_crit200 = np.fromfile(f, dtype=np.float32, count=ngroups)
    self.group_m_tophat200 = np.fromfile(f, dtype=np.float32, count=ngroups)
    self.group_veldisp = np.fromfile(f, dtype=np.float32, count=ngroups)

    if have_veldisp:
      self.group_veldisp_mean200 = np.fromfile(f, dtype=np.float32, count=ngroups)
      self.group_veldisp_crit200 = np.fromfile(f, dtype=np.float32, count=ngroups)
      self.group_veldisp_tophat200 = np.fromfile(f, dtype=np.float32, count=ngroups)

    self.group_nsubs = np.fromfile(f, dtype=np.uint32, count=ngroups)
    self.group_firstsub = np.fromfile(f, dtype=np.uint32, count=ngroups)

    if nsubs > 0:
      self.sub_len = np.fromfile(f, dtype=np.uint32, count=nsubs)
      self.sub_offset = np.fromfile(f, dtype=np.uint32, count=nsubs)
      self.sub_grnr = np.fromfile(f, dtype=np.uint64, count=nsubs)
      self.sub_nr = np.fromfile(f, dtype=np.uint64, count=nsubs)
      self.sub_pos = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=nsubs)
      self.sub_vel = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=nsubs)
      self.sub_cm = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=nsubs)
      # group_cm  (N,3) - Center of mass of the subgroup
      # group_vel (N,3) - Velocity of the subgroup 
      # group_pos (N,3) - Spatial position within the periodic box 
      self.sub_spin = np.fromfile(f, dtype=np.dtype((np.float32,3)), count=nsubs)
      self.sub_veldisp = np.fromfile(f, dtype=np.float32, count=nsubs)
      # sub_veldisp N   - One-dimensional velocity dispersion of all the member particles/cells (the 3D dispersion divided by sqrt(3) ).
      self.sub_vmax = np.fromfile(f, dtype=np.float32, count=nsubs)
      self.sub_vmaxrad = np.fromfile(f, dtype=np.float32, count=nsubs)
      self.sub_halfmassrad = np.fromfile(f, dtype=np.float32, count=nsubs)
      #self.sub_shape = np.fromfile(f, dtype=np.dtype((np.float32,6)), count=nsubs)
      self.sub_ebind = np.fromfile(f, dtype=np.float32, count=nsubs)
      self.sub_pot = np.fromfile(f, dtype=np.float32, count=nsubs)
      #self.sub_profile = np.fromfile(f, dtype=np.dtype((np.float32,9)), count=nsubs)
      self.sub_parent = np.fromfile(f, dtype=np.uint32, count=nsubs)
      self.sub_idbm = np.fromfile(f, dtype=np.uint64, count=nsubs)

      curpos = f.tell()
      f.seek(0,os.SEEK_END)
      if curpos != f.tell():
        print("Warning: finished reading before EOF for file",filenum)
      f.close()  

    
def read_groups_subhalos(filename, blocks):
    dmpm = 0.03722953 #dark matter particle mass
    DATAALL = {}
    DATATYPE = {}
    for block in blocks: 
        DATAALL[block] = []
        DATATYPE[block] = []
    grp_= subfind_catalog(filename) 
    nh_ = np.shape(grp_.group_nr)[0] #number of groups
    ns_ = np.shape(grp_.sub_nr)[0] #number of subhalos 
    #group information
    for block in blocks:
        if block == 'group_mass':
                DATAALL[block]  = np.log10(grp_.group_len*dmpm)+10
                DATATYPE[block] = 'float32'
        elif block == 'sub_mass':
                DATAALL[block]  = np.log10(grp_.sub_len*dmpm)+10 
                DATATYPE[block] = 'float32'
        else:
                DATAALL[block]  = getattr(grp_, block) 
                DATATYPE[block] = getattr(grp_, block).dtype
    return nh_, ns_, DATAALL, DATATYPE

def read_groups_subhalos_split(filenames, blocks):
    '''
    The worker of read_groups_subhalos for parallel 
    '''
    DATA = {}; DATATYPE = {} 
    for block in blocks:  DATA[block] = []
    
    Ns = []; Nh = [] 
    for filename in filenames:
        Nh_, Ns_, DATA_, DATATYPE_ = read_groups_subhalos(filename, blocks) 
        for block in blocks: DATA[block].append(DATA_[block])
        Nh.extend([Nh_])
        Ns.extend([Ns_])
    Nh = np.sum(Nh, dtype = np.int64)
    Ns = np.sum(Ns, dtype = np.int64)
    DATATYPE = DATATYPE_
    blockin3D = ['group_cm', 'group_vel', 'group_pos', 'sub_pos', 'sub_vel', 'sub_cm', 'sub_spin']
    for block in blocks: 
        if block in blockin3D: 
            DATA[block] = np.vstack(DATA[block])
        else: 
            DATA[block] = np.hstack(DATA[block])
    return Nh, Ns, DATA, DATATYPE

def read_groups(snapnum, blocks, basedir_groups = '/home/cossim/Jiutian/M1000/groups/'): 
    '''
    extract halo/subhalo informations from all subfind group files with multiprocessing 

    Wirtten by Qingyang Li, Yizhou Gu

    Parameters: 
    ----------
    snapnum: int
          the number of snapshot
    blocks:  list
          the name of blocks
    basedir_groups: str 
          file path for subfind data (FoF)

    Returns:
    -------
    DATAALL: dict  
          datasets of blocks. For example, group_pos = DATAALL['group_pos'][:]
    DATATYPE: dict 
          data type of blocks
    
    Note also: 
    ---------
    Host halo and subhalo information can be obtained together.
    '''
        
    #initialize data
    DATAALL = {}
    for block in blocks: DATAALL[block] = []
    filedir   = basedir_groups + '/groups_'+  str(snapnum).zfill(3) + \
                '/subhalo_tab_'+ str(snapnum).zfill(3)  + '.'
    filenames = glob.glob(filedir + '*') 
    division  = len(filenames) # the divided  number
    filenames = [filedir + '%s'%d for d in range(division) ] 

    size        = mp.cpu_count() 
    
    if division == 0:
        print('No found in ' + filedir)
        exit()
    else:
        print('Reading %s groups divisions of %sth snapshot with %s threads'%(division, snapnum, size)) 

    #split task & start parallel
    task_splits = np.array_split(filenames, size)
    pool = mp.Pool(size) 
    res  = []
    for ii in range(size):
        r = pool.apply_async(read_groups_subhalos_split, args = (task_splits[ii], blocks) ) 
        res.append(r)
    pool.close() 
    pool.join()

    NSUB = []
    for r in res:
        Nh_, Ns_, DATA, DATATYPE = r.get() 
        for block in blocks: DATAALL[block].append(DATA[block])
        
    blockin3D = ['group_cm', 'group_vel', 'group_pos', 'sub_pos', 'sub_vel', 'sub_cm', 'sub_spin']
    for block in blocks: 
        if block in blockin3D: 
            DATAALL[block] = np.vstack(DATAALL[block])
        else: 
            DATAALL[block] = np.hstack(DATAALL[block])
    return DATAALL, DATATYPE
    
def collect4csstmock(snapnum, basedir_hbt, basedir_groups):
    '''
    Collect input data for construction of csstmock from Jiutian using multiprocessing. 
    Using 72 threads, it takes ~ 20 min (10 min for readin and 10 min for finding next halo id).

    Wirtten by Qingyang Li, Yizhou Gu

    Parameters:
    ----------
    snapnum: int 
          the number of snapshot
    basedir_hbt: str
          file path for hbt data
    basedir_groups: str
          file path for subfind data (FoF)

    Return:
    -------
    Nsubs: int 
           number of subhalos
    INPUT_FORT_LONG: ndarray 
           include all int64  type data information
    INPUT_FORT_DOUBLE: ndarray 
           include all float64 type data information

    Details of subhalo information:
    --------

    >>> INPUT_FORT_LONG (np.int64):
    >>> 1. 'host_id':   host halo ID in present snapshot
    >>> 2. 'rank':      subhalo order (rank==1: central halo (most massive); rank ==0: subhalos) 
    >>> 3. 'sub_id':    trackID 
    >>> 4. 'host_nextid': host halo ID in next snapshot (-99 means not found)

    >>> INPUT_FORT_DOUBLE (np.float64): 
    >>> 1-3. 'sub_pos': position (x, y, z) 
    >>> 4-6. 'sub_velocity': velocity (vx, vy, vz): 
    >>> 7.   'host_mass': host halo mass
    >>> 8.   'sub_mass':  subhalo mass:  
    >>> 9.   'Peakmass': Peakmass
    >>> 10.  'sub_velmax': maximum circular velocity of subhalo
    >>> 
    >>> Note: 
    >>> mass in unit of Msun/h using log values 
    >>> cordinates in unit of Mpc/h 
    >>> velocity in unit of km/s
    '''

    start_time = time.time() 
#
#--- reading the present snapshot ---
#
    # subhalo information from HBT 
    blocks = ['ComovingAveragePosition', 'PhysicalAverageVelocity', 'LastMaxMass','VmaxPhysical','Mbound', 'Nbound', 'Rank', 'TrackId', 'HostHaloId']
    DATACOLLECT_hbt, DATATYPE_hbt = read_hbt(snapnum, blocks, basedir_hbt)
    # host halo information from subfind
    blocks = ['group_nr', 'group_mass']
    DATACOLLECT_groups, DATATYPE_groups   = read_groups(snapnum, blocks, basedir_groups) 

#
#--- reading next snapshot ---
#
    if snapnum != 127: 
        #subhalo information from HBT
        blocks = ['TrackId', 'HostHaloId']
        DATACOLLECT_hbt_next, DATATYPE_hbt_next = read_hbt(snapnum+1, blocks, basedir_hbt)
       
        #host halo information from subfind
        blocks = ['group_nr'] 
        DATACOLLECT_groups_next, DATATYPE_groups_next = read_groups(snapnum+1, blocks, basedir_groups)
    end_time = time.time()
    
    if snapnum != 127: 
       print('reading %sth and %sth snapshots time is %.3f s' %(snapnum, snapnum+1, end_time - start_time))
    if snapnum == 127: 
       print('reading %sth snapshots time is %.3f s' %(snapnum, end_time - start_time))
#---- finish reading, ~ 10 min.  

#==============================================================================
#
#--- prepare the data array for the construction of CSST mock 
#
    # the number of subhalo from HBT 
    Nsubs    = DATACOLLECT_hbt['TrackId'][:].shape[0] 
    arg_hbt  = np.argsort(DATACOLLECT_hbt['TrackId'][:])
    if snapnum != 127:
        arg_hbt_next = np.argsort(DATACOLLECT_hbt_next['TrackId'][:])
    
    
    # Data array which is Long type

    INPUT_FORT_LONG   = np.empty((Nsubs, 4),  dtype = np.int64)
    INPUT_FORT_LONG[:,0] = DATACOLLECT_hbt['TrackId'][arg_hbt]
    INPUT_FORT_LONG[:,1] = DATACOLLECT_hbt['Rank'][arg_hbt]
    idx_host   = DATACOLLECT_hbt['HostHaloId'][:] # hbt gives the index of host halos log Mh [Msub/h]
    hosthaloid = DATACOLLECT_groups['group_nr'][idx_host]
    idx_field = np.where(idx_host == -1)[0]
    hosthaloid[idx_field] = -1 # for field subhalos
    INPUT_FORT_LONG[:,2] = hosthaloid[arg_hbt]
    
    # determine the host halo id in next snapshot
    if snapnum != 127: 
        idx_host_next  = DATACOLLECT_hbt_next['HostHaloId'][:]
        next_groups_id = DATACOLLECT_groups_next['group_nr'][idx_host_next]
        idx_field_next = np.where(idx_host_next == -1)[0]
        next_groups_id[idx_field_next] = -1 # for field subhalos
        
        nexthaloid_time = time.time()
        INPUT_FORT_LONG[:,3] = nexthaloid(  DATACOLLECT_hbt['TrackId'][arg_hbt], \
                                            DATACOLLECT_hbt_next['TrackId'][arg_hbt_next], \
                                            next_groups_id[arg_hbt_next] )
    if snapnum == 127: INPUT_FORT_LONG[:,3] = -1
    
    # print('time for finding next host halo id is %.3f s' %(time.time() - nexthaloid_time) )

    # Data array which is Double type

    INPUT_FORT_DOUBLE        = np.empty((Nsubs, 10), dtype = np.float64)
    INPUT_FORT_DOUBLE[:,0:3] = DATACOLLECT_hbt['ComovingAveragePosition'][arg_hbt]
    INPUT_FORT_DOUBLE[:,3:6] = DATACOLLECT_hbt['PhysicalAverageVelocity'][arg_hbt]
    hosthalomass   = DATACOLLECT_groups['group_mass'][idx_host]
    hosthalomass[idx_field] = -1
    INPUT_FORT_DOUBLE[:,6]   = hosthalomass[arg_hbt]
    INPUT_FORT_DOUBLE[:,7]   = np.log10(DATACOLLECT_hbt['Mbound'][arg_hbt]) + 10
    INPUT_FORT_DOUBLE[:,8]   = np.log10(DATACOLLECT_hbt['LastMaxMass'][arg_hbt]) + 10 
    INPUT_FORT_DOUBLE[:,9]   = DATACOLLECT_hbt['VmaxPhysical'][arg_hbt] 
    
    end_time = time.time()
    print('reading the snapshot takes %.3f s in total' %(end_time - start_time) )
    return Nsubs, INPUT_FORT_LONG, INPUT_FORT_DOUBLE


def readsnapshot(snapnum, basedir_hbt, basedir_groups):
    '''
    Collect input data for construction of csstmock from Jiutian using multiprocessing. 
    Using 72 threads, it takes ~ 20 min (10 min for readin and 10 min for finding next halo id).

    Wirtten by Qingyang Li, Yizhou Gu

    Parameters:
    ----------
    snapnum: int 
          the number of snapshot
    basedir_hbt: str
          file path for hbt data
    basedir_groups: str
          file path for subfind data (FoF)

    Return:
    -------
    Nsubs: int 
           number of subhalos
    INPUT_FORT_LONG: ndarray 
           include all int64  type data information
    INPUT_FORT_DOUBLE: ndarray 
           include all float64 type data information

    Details of subhalo information:
    --------

    >>> INPUT_FORT_LONG (np.int64):
    >>> 1. 'host_id':   host halo ID in present snapshot
    >>> 2. 'rank':      subhalo order (rank==1: central halo (most massive); rank ==0: subhalos) 
    >>> 3. 'sub_id':    trackID 
    >>> 4. 'host_nextid': host halo ID in next snapshot (-99 means not found)

    >>> INPUT_FORT_DOUBLE (np.float64): 
    >>> 1-3. 'sub_pos': position (x, y, z) 
    >>> 4-6. 'sub_velocity': velocity (vx, vy, vz): 
    >>> 7.   'host_mass': host halo mass
    >>> 8.   'sub_mass':  subhalo mass:  
    >>> 9.   'Peakmass': Peakmass
    >>> 10.  'sub_velmax': maximum circular velocity of subhalo
    >>> 
    >>> Note: 
    >>> mass in unit of Msun/h using log values 
    >>> cordinates in unit of Mpc/h 
    >>> velocity in unit of km/s
    '''

    start_time = time.time() 
#
#--- reading the present snapshot ---
#
    # subhalo information from HBT 
    blocks = ['ComovingAveragePosition', 'PhysicalAverageVelocity', 'LastMaxMass','VmaxPhysical','Mbound', 'Nbound', 'Rank', 'TrackId', 'HostHaloId']
    DATACOLLECT_hbt, DATATYPE_hbt = read_hbt(snapnum, blocks, basedir_hbt)
    # host halo information from subfind
    blocks = ['group_nr', 'group_mass']
    DATACOLLECT_groups, DATATYPE_groups   = read_groups(snapnum, blocks, basedir_groups) 
    end_time = time.time() 

    print('reading %sth snapshots time is %.3f s' %(snapnum, end_time - start_time))
#---- finish reading, ~ 10 min.  

#==============================================================================
#
#--- prepare the data array for the construction of CSST mock 
#
    # the number of subhalo from HBT 
    Nsubs    = DATACOLLECT_hbt['TrackId'][:].shape[0] 
    #arg_hbt  = np.argsort(DATACOLLECT_hbt['TrackId'][:])

    # Data array which is Long type

    INPUT_FORT_LONG   = np.empty((Nsubs, 3),  dtype = np.int64)
    idx_host   = DATACOLLECT_hbt['HostHaloId'][:] # hbt gives the index of host halos log Mh [Msub/h]
    hosthaloid = DATACOLLECT_groups['group_nr'][idx_host]
    idx_field = np.where(idx_host == -1)[0]
    hosthaloid[idx_field] = -1                    # for field subhalos  
    INPUT_FORT_LONG[:,0] = hosthaloid[:]
    INPUT_FORT_LONG[:,1] = DATACOLLECT_hbt['Rank'][:]
    INPUT_FORT_LONG[:,2] = DATACOLLECT_hbt['TrackId'][:]

    # Data array which is Double type

    INPUT_FORT_DOUBLE        = np.empty((Nsubs, 10), dtype = np.float64)
    INPUT_FORT_DOUBLE[:,0:3] = DATACOLLECT_hbt['ComovingAveragePosition'][:]
    INPUT_FORT_DOUBLE[:,3:6] = DATACOLLECT_hbt['PhysicalAverageVelocity'][:]
    hosthalomass   = DATACOLLECT_groups['group_mass'][idx_host]
    hosthalomass[idx_field]  = -1
    INPUT_FORT_DOUBLE[:,6]   = hosthalomass[:]
    INPUT_FORT_DOUBLE[:,7]   = np.log10(DATACOLLECT_hbt['Mbound'][:]) + 10
    INPUT_FORT_DOUBLE[:,8]   = np.log10(DATACOLLECT_hbt['LastMaxMass'][:]) + 10 
    INPUT_FORT_DOUBLE[:,9]   = DATACOLLECT_hbt['VmaxPhysical'][:]
    end_time = time.time()
    print('reading the snapshot takes %.3f s in total' %(end_time - start_time) )
    return Nsubs, INPUT_FORT_LONG, INPUT_FORT_DOUBLE



def nexthaloid(TrackId, TrackIdnext, GroupsIdnext): 
    '''
    determine the host halo id in next snapshot. 
    
    Parameters:
    ----------
    TrackId: ndarray, int64 
            TrackId in present snapshot
    TrackIdnext: ndarray, int64
            TrackId in next snapshot
    GroupsIdnext: ndarray, int64
            GroupsId (host halo id) in next snapshot from subfind data (FoF)
    
    Return:
    -------
    GroupsId: ndarray, int64
            the corresponding GroupsId (host halo id) in next snapshot or -99 (no match) 
    '''
    
    # finding subhalo id in next snapshot
    idx_in_next    = np.in1d( TrackIdnext, TrackId, assume_unique = True)
    TrackId_match  = TrackIdnext[idx_in_next]  # Using TrackId to match
    GroupsID_match = GroupsIdnext[idx_in_next]
    
    # finding host halo id in next snapshot correspond to the present subhalos
    ind_in_prev    = np.in1d(  TrackId, TrackId_match, assume_unique = True) 
    GroupsId       = np.zeros( TrackId.shape[0], dtype = np.int64) - 99 
    GroupsId[ np.nonzero(ind_in_prev)[0] ] = GroupsID_match
    return GroupsId



import pandas as pd
import numpy as np
from pathlib import Path

def read_edition0(props, save=False): 
    '''
    Read the light-cone of Jiu-tian generated by Zhenlin Tan. 
    -->  a small, test edition of Jiutian lightcone.
    parameters:
    -----------
    
    
    return: 
    -------
    
    This script is written by Zhenlin Tan
    '''
    t1 = time.time() 

    output_list = ['trackID', 'hosthaloID', 'rank', 'posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'v_lfs',
                   'shMbound', 'd_comoving', 'ra', 'dec', 'vmax', 'PeakMass', 'PeakVmax', 'Mvir','Rvir',
                   'shBoundM200Crit','M200mean', 'redshift_true', 'redshift_obs', 'shMbound_at_ac', 'snapNum']
    global path, filepath, a, ret, names
    path = '/home/zltan/9tian-lightcone/testedition/'
    galbins = 2689
    print('Total files will be read: ', galbins)
    if props[0] == 'all':
        names = output_list
    else:
        names = props
    print('Properties will be read: ', names)

    for i in range(len(names)):  # 创建变量
        if names[i] in output_list:
            exec('%s = []' % names[i], globals())
        else:
            print('Input property %s is not in prop_list !' % names[i])
            return 0

#------------------------------------------------------------------------------

    for i in range(0, galbins):
        filepath = path + '9tmock_z05_%d.hdf5' % i
        if Path(filepath).is_file(): 
            a = pd.read_hdf(filepath) 
            # print( a )
            if i == 0: 
                for j in range(len(names)):
                    exec("%s = np.array(a['%s'])" % (names[j], names[j]) , globals() )
            elif len(a): 
                for j in range(len(names)):
                    exec("%s = np.array(np.concatenate((%s, a['%s']), axis=0))" % (names[j], names[j], names[j]), globals())
        else:
            print("File '%s' not exit !" % filepath)
        if i % 500 == 0 and len(names) and i != 0:
            print('Now %d files are read.' % i)
            exec("print(i, '/', galbins, '\t length now = ', len(%s))" % names[0])
    if save:
        for i in range(len(names)):
            exec("np.save('9tlc0_%s.npy', %s)" % (names[i], names[i]))  ###
            exec('''print('Save "9tlc0_%s.npy" done, length = ', len(%s))''' % (names[i], names[i])) 

    test = ','.join(str(i) for i in names)
    exec('ret = pd.DataFrame(np.array([%s]).T, columns=names)' % (test), globals())
    print('Data shaped : ', ret.shape)
    print('Read data done !')
    t2 = time.time() 
    print('It takes %s s.'% (t2 - t1) ) 
    return ret 

def lightcone4csstmock():
    
    ''' 
    Read the Jiutian lightcone for CSSTMOCK. A small, test edition of Jiutian lightcone by Zhenlin Tan.
    It takes ~600s to readin. 
    
    Provided by Zhenlin Tan, Jiaxin Han

    Return:
    -------
    Nsubs: int 
           number of subhalos
    INPUT_FORT_LONG: ndarray 
           include all int64   type data information
    INPUT_FORT_DOUBLE: ndarray 
           include all float64 type data information

    Details of lightcone data:
    --------
    
    >>> INPUT_FORT_LONG (np.int64):
    >>> 1. 'host_id':   host halo ID, from fof, 'group_nr'
    >>> 2. 'rank':      subhalo order
    >>> 3. 'sub_id':    'TrackID' 
    
    >>> INPUT_FORT_DOUBLE (np.float64): 
    >>> 1.   'ra': Right Ascension
    >>> 2.   'dec': Declination 
    >>> 3.   'd_comoving': comoving distance
    >>> 4.   'redshift_true': cosmology redshift, z_cos
    >>> 5.   'redshift_obs': observed redshift, z_obs 
    >>> 6.   'v_lfs': line-of-sight velocity 
    >>> 7.   'group_mass':  host halo mass from fof, 'group_mass'
    >>> 8.   'shMbound': subhalo mass 
    >>> 9.   'PeakMass': peak mass  
    >>> 10.   'vmax':     max circular velocity 
    >>> 
    >>> Note: 
    >>> mass in unit of Msun/h using log values 
    >>> cordinates in unit of Mpc/h 
    >>> velocity in unit of km/s
    
    '''
# prop_list = ['trackID', 'hosthaloID', 'rank', 'posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'v_lfs',
#              'shMbound', 'd_comoving', 'ra', 'dec', 'vmax', 'PeakMass', 'PeakVmax', 'Mvir','Rvir',
#              'shBoundM200Crit','M200mean', 'redshift_true', 'redshift_obs', 'shMbound_at_ac', 'snapNum']

    props = ['trackID', 'snapNum', 'hosthaloID', 'rank', \
             'ra', 'dec', 'd_comoving', \
             'redshift_true', 'redshift_obs', 'v_lfs', \
             'Mvir', 'shMbound', 'PeakMass', 'vmax'] 
    ret = read_edition0(props, save=False)
    ret = ret.sort_values(by=['trackID', 'snapNum'])

    group_path = '/home/zltan/9tian-lightcone/testedition/group_patch.hdf5'
    group_props = pd.read_hdf(group_path)
#    print(ret.iloc[0:10])
#    print(group_props.iloc[0:10]) 

    Nsubs             = ret['trackID'].shape[0]  
    # print(ret['trackID'].shape )
    INPUT_FORT_LONG   = np.empty((Nsubs, 3),  dtype = np.int64)
    INPUT_FORT_LONG[:,0]   = group_props['group_nr'][:Nsubs] 
    INPUT_FORT_LONG[:,1]   = ret['rank']
    INPUT_FORT_LONG[:,2]   = ret['trackID'] 

    # indx = 7324 == INPUT_FORT_LONG[:,0]
    # print(INPUT_FORT_LONG[indx,1]) 


    INPUT_FORT_DOUBLE  = np.empty((Nsubs, 10),  dtype = np.float64)
    INPUT_FORT_DOUBLE[:,0]  = ret['ra']
    INPUT_FORT_DOUBLE[:,1]  = ret['dec']
    INPUT_FORT_DOUBLE[:,2]  = ret['d_comoving']
    INPUT_FORT_DOUBLE[:,3]  = ret['redshift_true']
    INPUT_FORT_DOUBLE[:,4]  = ret['redshift_obs']
    INPUT_FORT_DOUBLE[:,5]  = ret['v_lfs']
    INPUT_FORT_DOUBLE[:,6]  = np.array(group_props['group_mass'])[:Nsubs] 
    INPUT_FORT_DOUBLE[:,7]  = np.log10(np.array(ret['shMbound']).astype('float64') ) + 10 
    INPUT_FORT_DOUBLE[:,8]  = np.log10(np.array(ret['PeakMass']).astype('float64') ) + 10 
    INPUT_FORT_DOUBLE[:,9]  = ret['vmax']
    

    return Nsubs, INPUT_FORT_LONG,  INPUT_FORT_DOUBLE

if __name__ == '__main__': 
    #Nsubs, INPUT_FORT_LONG,  INPUT_FORT_DOUBLE = lightcone4csstmock() 
    basedir = '/home/cossim/Jiutian/M1000/' 
    basedir_hbt    = basedir + '/hbt/'
    basedir_groups = basedir + '/groups/'
    Nsubs, INPUT_FORT_LONG,  INPUT_FORT_DOUBLE = collect4csstmock(127, basedir_groups=basedir_groups, basedir_hbt= basedir_hbt) 
    np.save('INPUT_FORT_LONG', INPUT_FORT_LONG) 
    np.save('INPUT_FORT_DOUBLE', INPUT_FORT_DOUBLE) 
    
if __name__ == '__main__0': 
    snapnum = 127 
    basedir = '/home/cossim/Jiutian/M1000/' 
    basedir_hbt    = basedir + '/hbt/'
    basedir_groups = basedir + '/groups/'
    start_time = time.time()
    Nsubs, INPUT_FORT_LONG, INPUT_FORT_DOUBLE = collect4csstmock(snapnum, basedir_hbt, basedir_groups) 
    end_time = time.time()
    print('The time of collect4csstmock for %s snapshot is %.3f s' %(snapnum, end_time - start_time))
    print(Nsubs) 
    print(INPUT_FORT_LONG.shape) 
    print(INPUT_FORT_DOUBLE.shape)
    np.save('INPUT_FORT_LONG', INPUT_FORT_LONG) 
    np.save('INPUT_FORT_DOUBLE', INPUT_FORT_DOUBLE) 
