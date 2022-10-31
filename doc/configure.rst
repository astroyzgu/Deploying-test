pipeline configure 
################## 

.. code-block:: python  

    # csstmock.init 
    [basic]
    ncore   = 72   
    cosmo
    [lightcone]
    simudir = '/home/cossim/Jiutian/M1000/' 
    simuname= 'jiutian' # ELUCID
    subhalo = 'hbt'  # ['hbt' or 'fof']  
    # input:  
    # x,y,z,vx,vy,vz, Nin, hid, subhid, next_hid
    # rsmin, rsmax
    # boxsize = 1000 # Mpc/h 
    # output:
    # ra, dec, r, z_cos, z_obs, vr  
    # hid, subhid, next_hid
    # Nout 
    lumfunc = './data/lumfunc_example.txt' # select the luminosity function 
    # output: logL, logn 
    [mockgalaxy] 
    skip = 'False'
    [analysis]
    skip = 'False'
    [outputs]
    outputs = './'

