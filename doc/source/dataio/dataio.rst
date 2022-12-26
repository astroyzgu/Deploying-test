.. currentmodule:: csstmock.dataio.io_jiutian

:mod:`dataio` module 
========================

:mod:`io_jiutian` -- dataio of the jiutian simulation 
----------------------------------

The Jiutian Main Simulation is the high-resolution N-body simulations of strand $\rm \Lambda CDM$ cosmology with $\Omega_m=0.3111$, $\Omega_\Lambda=0.6889$. 
This simulation box, running with Gadget3, is 1 Gpc/h with $6144^3$ dark matter particles of $3.72 \times 10^8 M_{\odot}/h$.  

- :func:`read_hbt` read subhalo information from all :doc:`jiutian_hbt` subhalo files.  
- :func:`read_groups` extract halo/subhalo informations from all :doc:`jiutian_subfind`  group files.

Example 
^^^^^^^^^
.. code-block:: python  

    from csstmock.dataio.io_jiutian import read_hbt, read_groups 
    snapnum = 127 
    basedir = '/home/cossim/Jiutian/M1000/' 
    basedir_hbt    = basedir + '/hbt/'
    basedir_groups = basedir + '/groups/' 
    blocks = ['ComovingAveragePosition', 'PhysicalAverageVelocity', \
              'LastMaxMass','VmaxPhysical', \ 
              'Mbound', 'Nbound', 'Rank', \
              'TrackId', 'HostHaloId']
    DATACOLLECT_hbt, DATATYPE_hbt = read_hbt(snapnum, blocks, basedir_hbt)

    # host halo information from subfind
    blocks = ['group_nr', 'group_mass']
    DATACOLLECT_groups, DATATYPE_groups = read_groups(snapnum, blocks, basedir_groups)     
