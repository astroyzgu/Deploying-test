import os 
import numpy as np 
import healpy as hp 
def skycov_healpy(survey): 
    LOCATION   = os.path.dirname(os.path.abspath(__file__)) 
    if survey == 'desidr9': 
        wmap_healpy = os.path.join(LOCATION, '../../', './database_builder/skycov/skycov256-desidr9.npy')
        w        = np.load(wmap_healpy)
        w        = w.astype('float64') 
        nside    = np.sqrt( np.shape(w)[0]/12.0 )
        if not nside in 2**np.arange(0, 30): 
            raise ValueError( "Running as HEALPY; %s pixels are found."% indx.shape[0] + \
            "However, corresponding nside (%s) must equal to 2^N, N = 1,2,3,4,...,29. \n"% nside + \
            "Thus, corresponding nside must larger than nside >= 2^%s (%s)"%( int( np.log2(nside) ), 2**int( np.log2(nside) ) ) ) 
        nside = int(nside)
    if (isinstance(survey, int)): 
        nside = survey 
        if nside in 2**np.arange(0, 30): 
            w = np.arange(12*nside*nside) 
        else: 
            raise ValueError("nside (%s) must equal to 2^N, N = 1,2,3,4,...,29. \n"% nside)

    return w, nside

def assignwht_healpy(x, y, w): 
    w        = np.array(w)
    nside    = np.sqrt( np.shape(w)[0]/12.0 )
    if not nside in 2**np.arange(0, 30): 
            raise ValueError( "Running as HEALPY; %s pixels are found."% indx.shape[0] + \
            "However, corresponding nside (%s) must equal to 2^N, N = 1,2,3,4,...,29. \n"% nside + \
            "Thus, corresponding nside must larger than nside >= 2^%s (%s)"%( int( np.log2(nside) ), 2**int( np.log2(nside) ) ) ) 
    nside = int(nside) 
    ipix = hp.ang2pix(nside, x, y, lonlat = True)
    return w[ipix] 

if '__main__' == __name__: 
    x = np.array([1,40])
    y = np.array([1,40]) 
    wht, nside = skycov_healpy('desidr9') 
    w          = assignwht_healpy(x,y, wht) 
    print(w)
    hp.mollview(wht, rot = [118, 0, 0])
    hp.projscatter(x, y, lonlat = True, )