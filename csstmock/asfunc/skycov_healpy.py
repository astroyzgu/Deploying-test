import os 
import numpy as np 
import healpy as hp 
from astropy.coordinates import SkyCoord

def skycov_healpy(survey): 
    LOCATION   = os.path.dirname(os.path.abspath(__file__)) 
    survey_available =   ['desidr9', 'hscdr3', 'csstv0']
    if not survey in survey_available: 
        raise ValueError('%s is not available. Only the following is avaible:', survey_available )

    if survey == 'desidr9': 
        wmap_healpy = os.path.join(LOCATION, '../../', './database_builder/skycov/skycov256-desidr9.npy')
        
    if survey == 'hscdr3': 
        wmap_healpy = os.path.join(LOCATION, '../../', './database_builder/skycov/skycov512-hscdr3.npy') 

    if survey == 'csstv0': 
        wmap_healpy = os.path.join(LOCATION, '../../', './database_builder/skycov/skycov512-csstv0.npy') 
    
    pix      = np.load(wmap_healpy).astype('int64')
    nside    = wmap_healpy.split('/')[-1].split('-')[0].replace('skycov', '')
    nside    = int(nside)
    if not nside in 2**np.arange(0, 30):
            raise ValueError( "Running as HEALPY; %s pixels are found."% indx.shape[0] + \
            "However, corresponding nside (%s) must equal to 2^N, N = 1,2,3,4,...,29. \n"% nside + \
            "Thus, corresponding nside must larger than nside >= 2^%s (%s)"%( int( np.log2(nside) ), 2**int( np.log2(nside) ) ) ) 
            
    w = np.arange(12*nside*nside).astype('float64')*0.0
    w[pix] = 1.0
    
    npix = len(pix)
    pixarea = hp.nside2pixarea(nside, degrees= True) 
    print( 'sky coverage is %.2f deg^2 using nside = %s (pixarea = %9.3e deg^2/pixel)'%(pixarea*npix, nside, pixarea) )
    ra, dec = hp.pix2ang(nside, pix, lonlat = True)
    coord   = SkyCoord(ra, dec, unit = 'deg', frame='icrs')  # using degrees directly
    npix = int( np.sum( coord.galactic.b < 0) )
    print( 'South (b<0) coverage = %.2f deg^2'%(pixarea*npix ) )
    npix = int( np.sum( coord.galactic.b > 0) )
    print( 'North (b>0) coverage = %.2f deg^2'%(pixarea*npix ) )
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
import matplotlib.pyplot as plt

if '__main__' == __name__: 
    x = np.array([1,40])
    y = np.array([1,40]) 
    #wht, nside = skycov_healpy('desidr9') 
    #wht, nside = skycov_healpy('hscdr3') 
    wht, nside = skycov_healpy('csstv0') 
    ipix = np.arange(12*nside*nside)
    ipix = ipix[wht==1]
    w          = assignwht_healpy(x,y, wht) 
    print(w)
    hp.mollview(wht, rot = [118, 0, 0])
    hp.projscatter(x, y, lonlat = True, )
    plt.show() 




