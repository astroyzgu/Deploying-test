Fortran Interface
==========

USAGE
----------

.. code-block:: shell  

    # 编译动态链接库
    module load compiler/intel-2018 # 加载ifort(gravity上)
    python build.py                 # 生成静态链接库 libcsstplugin.a 
    gfortran -fPIC -shared -o libcsstmock.so io_dummy.f io_jiutian.f libcsstplugin.a # 合并为动态链接库
    
    # 添加路径:
    export PYTHONPATH=$PYTHONPATH:$(pwd)
    export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(pwd)
    
    # 编译程序 readsnapshot_jiutian的例子
    ifort -o ./test1 test1.f -lcsstmock
    ./test1
    
    # 编译程序 readlightcone_jiutian的例子
    ifort -o ./test2 test2.f -lcsstmock
    ./test2


subroutine 
----------
In current stage, the following subroutine is embedded:

- readsnap_jiutian_, read one snapshot from jiutian simulation 
- skycov_, assign weight [0 or 1] of a given survey ['desidr9', 'hscdr3', 'csstv0']. 

readsnap_jiutian
^^^^^^^^^^
.. code-block:: Fortran 

   subroutine readsnap_jiutian(snap, darr, larr, nd1, ngalmax, ngal) 

   Input: 
   snap(= 0-127), the number of snapshot 
   nd1(= 10), the size of dimension 1 of larr
   ngalmax(= 500000000), the size of dimension 2 of larr and darr. 
           Namely, the max number is allowed. 

   Output: 
   larr, int*8, the size of 2d array is (3, ngalmax)
    1. 'host_id': host halo ID in present snapshot
    2. 'sub_id':  subhalo ID
    3. 'host_nextid': host halo ID in next snapshot (-99 means not found)

   darr, real*8, the size of 2d array is (10, ngalmax) 
    1-3. 'sub_pos': position (x, y, z)
    4-6. 'sub_velocity': velocity (vx, vy, vz):
    8. 'host_mass': host halo mass
    8. 'sub_mass': subhalo mass:
    9. 'PeakMass': peak mass
    10. 'sub_velmax': maximum circular velocity of subhalo

   ngal, the real number which are loaded, 
         e.g., 460267974 subhalos at 127th snaphot although ngalmax is 5E8. 

skycov
^^^^^^^^^^
.. code-block:: Fortran

    subroutine skycov(x, y, w, survey, Ngalmax, Ngal)

    Input: 
    x, array-like, R.A. 
    y, array-like, Dec.
    survey, Optional. ['desidr9', 'hscdr3', 'csstv0'] 
    Ngalmax, the size (max number) of x, y, w is allowed. 
    Ngal, assign weight only for the first Ngal galaxies.  

    Output: 
    w, if within the survey, then 1; else, 0.   
