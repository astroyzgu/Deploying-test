Installation and dependecies 
=========

Installation 
---------
INSTALL the package

.. code-block:: bash  

   # 1. Copy the package:

   cp cssstmock-dev ./to/your/path/
   cd cssstmock-dev

   # 2. python part:

   #--- If on Gravity, 
   #--- just load two modules before INSTALL:

   module purge
   conda deactivate
   module load anaconda
   module load compiler/intel-2018  # 加载ifort(gravity上)
   source activate python3          # 激活gravity上自带的pythob3环境

   #--- If not on Gravity, 
   #--- build the python environments before INSTALL (not now):

   conda create -n csstmock  python=3.9 healpy
   conda activate csstmock
   conda install mamba -c conda-forge
   mamba install cffi ipykernel h5py pytables
   mamba install astropy matplotlib pandas

   #--- Then, INSTALL csstmock pipeline of PYTHON package:

   python setup.py develop --user

   # 3. Fortran part: 

   #--- Install the dynnamic libaray for fortran (also see the libcsstmock/readme.md)

   cd libcsstmock   # 移动到动态链接库的目录
   python build.py  # 生成静态链接库 libcsstplugin.a
   gfortran -fPIC -shared -o libcsstmock.so io_dummy.f io_jiutian.f libcsstplugin.a # 合并为动态链接库

   #--- 添加当前路径入库:
   export PYTHONPATH=$PYTHONPATH:$(pwd)
   export LIBRARY_PATH=$LIBRARY_PATH:$(pwd)
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
   export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$(pwd)
 
