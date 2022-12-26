import multiprocessing as mp 
import numpy as np 
import ctypes 
import time
import glob 
def mp2np(mp_array, shape = None):
    '''Given a multiprocessing.Array, returns an ndarray pointing to
    the same data.'''
    # support SynchronizedArray:
    _ctypes_to_numpy = {
    ctypes.c_char   : np.dtype(np.uint8),
    ctypes.c_wchar  : np.dtype(np.int16),
    ctypes.c_byte   : np.dtype(np.int8),
    ctypes.c_ubyte  : np.dtype(np.uint8),
    ctypes.c_short  : np.dtype(np.int16),
    ctypes.c_ushort : np.dtype(np.uint16),
    ctypes.c_int    : np.dtype(np.int32),
    ctypes.c_uint   : np.dtype(np.uint32),
    ctypes.c_long   : np.dtype(np.int64),
    ctypes.c_ulong  : np.dtype(np.uint64),
    ctypes.c_float  : np.dtype(np.float32),
    ctypes.c_double : np.dtype(np.float64)}
    if not hasattr(mp_array, '_type_'):
        mp_array = mp_array.get_obj()
    dtype = _ctypes_to_numpy[mp_array._type_]
    result = np.frombuffer(mp_array, dtype)
    if shape is not None:
        result = result.reshape(shape)
    return np.asarray(result)

def np2mp(np_array, lock = False):
    '''Generate an 1D multiprocessing.Array containing the data from
    the passed ndarray.  The data will be *copied* into shared
    memory.'''
    _ctypes_to_numpy = {
    ctypes.c_char   : np.dtype(np.uint8),
    ctypes.c_wchar  : np.dtype(np.int16),
    ctypes.c_byte   : np.dtype(np.int8),
    ctypes.c_ubyte  : np.dtype(np.uint8),
    ctypes.c_short  : np.dtype(np.int16),
    ctypes.c_ushort : np.dtype(np.uint16),
    ctypes.c_int    : np.dtype(np.int32),
    ctypes.c_uint   : np.dtype(np.uint32),
    ctypes.c_long   : np.dtype(np.int64),
    ctypes.c_ulong  : np.dtype(np.uint64),
    ctypes.c_float  : np.dtype(np.float32),
    ctypes.c_double : np.dtype(np.float64)}
    _numpy_to_ctypes = dict(zip(_ctypes_to_numpy.values(), _ctypes_to_numpy.keys()))
    array1d = np_array.ravel(order = 'A')
    try:
        c_type = _numpy_to_ctypes[array1d.dtype]
    except KeyError:
        c_type = _numpy_to_ctypes[np.dtype(array1d.dtype)]
    mp_array   = mp.Array( c_type, array1d.size, lock = lock)
    mp2np(mp_array)[:] = array1d
    return mp_array

def read_groups_division(filename): 
    data = np.random.normal(0,1,size=(10,3))
    return data 

def read_groups_divisions(filenames):
    data = []
    for filename in filenames: 
        data_ = read_groups_division(filename)
        data.append(data_)
    data = np.vstack(data)
    return data

def read_groups(snapnum, basedir = '/home/cossim/Jiutian/M1000/groups/'):

    #------- check whether are files existed.
    
    filedir   = basedir + 'groups_'+ str(snapnum).zfill(3) \
                + '/subhalo_tab_'+ str(snapnum).zfill(3) + '.*'
    filenames = glob.glob( filedir )
    division  = len(filenames) # the divided  number
    if division == 0:
        print('No found in ' + filedir)
        exit()
    else:
        print('%s divisions are found in %sth snapshot'%(division, snapnum)) 

    #------- loading data 
    
    size        = mp.cpu_count() 
    filenames   = [ basedir + 'groups_'+ str(snapnum).zfill(3)    \
                + '/subhalo_tab_'+ str(snapnum).zfill(3) + str(d) \
                   for d in range(division)]  
    task_splits = np.array_split( filenames, size)
    print('Runing with %s threads'% size) 
    pool = mp.Pool(size) 
    res  = []
    for ii in range(size):
        r = pool.apply_async( read_groups_divisions, args = (task_splits[ii],) ) 
        res.append(r)
    pool.close() 
    pool.join()
    return np.vstack([r.get() for r in res]) 
    
    
def read_hbt():
    pass
    return 

if __name__ == '__main__':  
    '''
    
    ''' 
    nsub         = 10000 
    d_np_arr2d    = np.empty( (nsub, 3), dtype = 'float64') 
    l_np_arr2d    = np.empty( (nsub, 3), dtype = 'int64') 
    d_mp_arr2d = np2mp(d_np_arr2d); d_np_arr2d = mp2np(d_mp_arr2d, d_np_arr2d.shape)
    l_mp_arr2d = np2mp(l_np_arr2d); l_np_arr2d = mp2np(l_mp_arr2d, l_np_arr2d.shape)
    basedir  = '/home/cossim/Jiutian/M1000/'
    snapnum  = 127  # the snapshot number 
    info_grp = read_groups(snapnum, basedir = basedir + '/groups/')
    print(info_grp.shape, info_grp.dtype) 
    info_hbt = read_hbt(snapnum, basedir = basedir + '/groups/')
