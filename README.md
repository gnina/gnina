gnina (pronounced NEE-na) is a fork of smina, which is a fork of AutoDock Vina.

gnina is **not** recommended for production use (*yet*) in molecular modeling tasks.  However, it *is* suitable as a platform for researching structure-based deep learning approaches as described in [our paper](http://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00740).  

Help
====
**Reminder: gnina is not yet intended for production use**.  However, if you would like to evaluate it or use it as a research platform, please [subscribe to our slack team](http://bits.csb.pitt.edu/slack). 

Installation
============

To install (Ubuntu 16.04):
```
apt-get install build-essential git wget libopenbabel-dev libboost-all-dev libeigen3-dev libgoogle-glog-dev libprotobuf-dev protobuf-compiler libhdf5-serial-dev libatlas-base-dev python-dev cmake librdkit-dev python-numpy
```

[Follow NVIDIA's instructions](http://docs.nvidia.com/cuda/cuda-installation-guide-linux/#axzz4TWipdwX1) to install the latest version of CUDA.  Or:

```
wget https://developer.nvidia.com/compute/cuda/8.0/Prod2/local_installers/cuda-repo-ubuntu1604-8-0-local-ga2_8.0.61-1_amd64-deb
dpkg -i cuda-repo-ubuntu1604-8-0-local_8.0.44-1_amd64-deb 
apt-get update
apt-get install cuda
```
```
git clone https://github.com/gnina/gnina.git
cd gnina
mkdir build
cd build
cmake ..
make
make install
```
# 

To install on (CentOS 7):

The program will not build in a computer with a gpu with computer capability < 3.5 unless
you setup force a different architecture.

_Add EPEL repository_
```
sudo yum  install epel-release
```
_[Follow NVIDIA's instructions](http://docs.nvidia.com/cuda/cuda-installation-guide-linux/#axzz4TWipdwX1) to install the latest version of CUDA.  Or:_
```
wget http://developer.download.nvidia.com/compute/cuda/repos/rhel7/x86_64/cuda-repo-rhel7-8.0.61-1.x86_64.rpm
sudo rpm -i cuda-repo-rhel7-8.0.61-1.x86_64.rpm
sudo yum clean all
sudo yum install cuda
```

_Install dependencies_
```
sudo yum  groupinstall 'Development Tools'
```
```
sudo yum install boost-devel.x86_64 eigen3-devel.noarch protobuf-compiler.x86_64 protobuf-devel.x86_64 hdf5-devel.x86_64 cmake git wget openbabel-devel.x86_64 openbabel.x86_64 leveldb-devel.x86_64 snappy-devel.x86_64 opencv-devel.x86_64 gflags-devel.x86_64 glog-devel.x86_64 lmdb-devel.x86_64 readline-devel.x86_64 zlib-devel.x86_64 bzip2-devel.x86_64 sqlite-devel.x86_64 python-devel.x86_64 numpy.x86_64 atlas-devel.x86_64 atlas.x86_64 atlas-static.x86_64
```
_Install cmake 3.8_  
The cmake installed by yum in CentOS 7 (cmake version 2.8.12.2) produce a lot of error. Is better if you use an updated version.
```
cd /home/$USER/bin
wget https://cmake.org/files/v3.8/cmake-3.8.0-Linux-x86_64.tar.gz
tar -xvf cmake-3.8.0-Linux-x86_64.tar.gz
export CMAKE_HOME=/home/$USER/bin/cmake-3.8.0-Linux-x86_64
export PATH=$CMAKE_HOME/bin:$PATH
```
_Install RDKit Release_2017_03_1 and compile gnina_  
_Install RDKit_  
Is better if we keep everything inside the gnina directory
```
cd /home/$USER/bin
git clone https://github.com/gnina/gnina.git
cd gnina
wget https://github.com/rdkit/rdkit/archive/Release_2017_03_1.tar.gz
tar -xvf Release_2017_03_1.tar.gz
cd Release_2017_03_1
export RDBASE=`pwd`
export LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH
mkdir build
cd build
```
**If you are using anaconda python the you need to check that all the python variables are set correctly or set them manually.**
```
export ANACONDA_PY_HOME=/home/$USER/anaconda
cmake -DPYTHON_EXECUTABLE=$ANACONDA_PY_HOME/bin/python -DPYTHON_INCLUDE_DIR=$ANACONDA_PY_HOME/include/python2.7 -DPYTHON_LIBRARY=$ANACONDA_PY_HOME/lib/libpython2.7.so -DPYTHON_NUMPY_INCLUDE_PATH=$ANACONDA_PY_HOME/lib/python2.7/site-packages/numpy/core/include ..
make
ctest
make install
```
**If you are using your CentOS python**
```
cmake ..
make 
ctest
make install
```
_Fix RDKit Libraries_  
Compiling RDKit will add the name of the package to the library.   
ex. libSmilesParse.so (UBUNTU Package) !=  libRDKitSmilesParse.so (Compiled in CentOS)  
We need to make additional links to resemble the UBUNTU names.
```
cd $RDBASE/lib
for i in $(ls -1 *.so.1.2017.03.1); do name=`basename $i .so.1.2017.03.1`; namef=`echo $name | sed 's/RDKit//g'`; ln -s $i ${namef}.so.1; ln -s ${namef}.so.1 ${namef}.so; done
```
_Continue with gnina compilation_  
We need to set the variable for the ATLAS libraries.  
Use lib**s**atlas.so for serial libraries or lib**t**atlas.so for threaded libraries.  

**If you are using anaconda python the you need to check that all the python variables are set correctly or set them manually.**
```
cd /home/$USER/bin/gnina
mkdir build
cd build
cmake -DPYTHON_EXECUTABLE=$ANACONDA_PY_HOME/bin/python -DPYTHON_INCLUDE_DIR=$ANACONDA_PY_HOME/include/python2.7 -DPYTHON_LIBRARY=$ANACONDA_PY_HOME/lib/libpython2.7.so -DAtlas_BLAS_LIBRARY=/usr/lib64/atlas/libtatlas.so -DAtlas_CBLAS_LIBRARY=/usr/lib64/atlas/libtatlas.so -DAtlas_LAPACK_LIBRARY=/usr/lib64/atlas/libtatlas.so ..
make 
make install
```
**If you are using your CentOS python**
```
cd /home/$USER/bin/gnina
mkdir build
cd build
cmake -DAtlas_BLAS_LIBRARY=/usr/lib64/atlas/libtatlas.so -DAtlas_CBLAS_LIBRARY=/usr/lib64/atlas/libtatlas.so -DAtlas_LAPACK_LIBRARY=/usr/lib64/atlas/libtatlas.so ..
make
make install
```

#
If you are building for systems with different GPUs, include `-DCUDA_ARCH_NAME=All`.  


Training
========

Scripts to aid in training new CNN models can be found at [https://github.com/gnina/scripts](https://github.com/gnina/scripts)
and sample models at [https://github.com/gnina/models](https://github.com/gnina/models).

The input layer should be a `MolGridData` layer.  For example:
```
layer {
  name: "data"
  type: "MolGridData"
  top: "data"
  top: "label"
  include {
    phase: TRAIN
  }
  molgrid_data_param {
    source: "TRAINFILE"
    batch_size:  20
    dimension: 23.5
    resolution: 0.5
    shuffle: true
    balanced: true
    random_rotation: true
    random_translate: 2
    root_folder: "/home/dkoes/CSAR/"
  }
}
```

This layer performs GPU-accelerated grid generation on-the-fly which means it can apply random rotations
and translations to the input (essential for training).  The input file (TRAINFILE) contains an example on each line, 
which consists of a label, a receptor file, and a ligand file:

```
1 set2/297/rec.gninatypes set2/297/docked_0.gninatypes # text after a hash is ignored
1 set2/297/rec.gninatypes set2/297/docked_1.gninatypes
1 set2/297/rec.gninatypes set2/297/docked_2.gninatypes 
1 set2/297/rec.gninatypes set2/297/docked_3.gninatypes 
0 set2/297/rec.gninatypes set2/297/docked_4.gninatypes 
0 set2/297/rec.gninatypes set2/297/docked_5.gninatypes 
...
```

Althoug the receptor and ligand can be specified as any normal molecular data file, we strongly recommend (for training at least)
that molecular structure files be converted to `gninatypes` files with the `gninatyper` executable.  These are much smaller files
that incur less I/O. Relative file paths will be prepended with the `root_folder` parameter in MolGridData, if applicable.

The provided models are templated with `TRAINFILE` and `TESTFILE` arguments, which the `train.py` script will substitue with 
actual files.  The `train.py` script can be called with a model and a prefix for testing and training files:

```
train.py -m models/refmodel3/refmodel3.model -p models/data/csar/all
```

This will perform cross-validation using the `alltrain[0-2].types` and `alltest[0-2].types` files.

There are quite a few options to `train.py` for modifying training:
```
usage: train.py [-h] -m MODEL -p PREFIX [-n NUMBER] [-i ITERATIONS] [-s SEED]
                [-t TEST_INTERVAL] [-o OUTPREFIX] [-g GPU] [-c CONT] [-k] [-r]
                [--avg_rotations] [--keep_best] [--dynamic] [--solver SOLVER]
                [--lr_policy LR_POLICY] [--step_reduce STEP_REDUCE]
                [--step_end STEP_END] [--step_when STEP_WHEN]
                [--base_lr BASE_LR] [--momentum MOMENTUM]
                [--weight_decay WEIGHT_DECAY] [--gamma GAMMA] [--power POWER]
                [--weights WEIGHTS]

Train neural net on .types data.

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL, --model MODEL
                        Model template. Must use TRAINFILE and TESTFILE
  -p PREFIX, --prefix PREFIX
                        Prefix for training/test files:
                        <prefix>[train|test][num].types
  -n NUMBER, --number NUMBER
                        Fold number to run, default is all
  -i ITERATIONS, --iterations ITERATIONS
                        Number of iterations to run,default 10,000
  -s SEED, --seed SEED  Random seed, default 42
  -t TEST_INTERVAL, --test_interval TEST_INTERVAL
                        How frequently to test (iterations), default 40
  -o OUTPREFIX, --outprefix OUTPREFIX
                        Prefix for output files, default <model>.<pid>
  -g GPU, --gpu GPU     Specify GPU to run on
  -c CONT, --cont CONT  Continue a previous simulation from the provided
                        iteration (snapshot must exist)
  -k, --keep            Don't delete prototxt files
  -r, --reduced         Use a reduced file for model evaluation if exists(<pre
                        fix>[_reducedtrain|_reducedtest][num].types)
  --avg_rotations       Use the average of the testfile's 24 rotations in its
                        evaluation results
  --keep_best           Store snapshots everytime test AUC improves
  --dynamic             Attempt to adjust the base_lr in response to training
                        progress
  --solver SOLVER       Solver type. Default is SGD
  --lr_policy LR_POLICY
                        Learning policy to use. Default is inv.
  --step_reduce STEP_REDUCE
                        Reduce the learning rate by this factor with dynamic
                        stepping, default 0.5
  --step_end STEP_END   Terminate training if learning rate gets below this
                        amount
  --step_when STEP_WHEN
                        Perform a dynamic step (reduce base_lr) when training
                        has not improved after this many test iterations,
                        default 10
  --base_lr BASE_LR     Initial learning rate, default 0.01
  --momentum MOMENTUM   Momentum parameters, default 0.9
  --weight_decay WEIGHT_DECAY
                        Weight decay, default 0.001
  --gamma GAMMA         Gamma, default 0.001
  --power POWER         Power, default 1
  --weights WEIGHTS     Set of weights to initialize the model with
```

User Grids
----------

In some cases it may be desirable to incorporate additional grid-based input
into the training data.  In this case it is necessary to pre-generate 
grids from the molecular data and user-supplied grids with `gninagrid` and use
the `NDimData` input layer.

```
layer {
  name: "data"
  type: "NDimData"
  top: "data"
  top: "label"
  include {
    phase: TRAIN
  }
  ndim_data_param {
    source: "TRAINFILE"
    batch_size: 10
    shape {
      dim: 34
      dim: 48
      dim: 48
      dim: 48
    }
    shuffle: true
    balanced: true
    rotate: 24
  }
}
```

Similar to the MolGrid layer, TRAINFILE contains an example on each line with a label and one or more binmap files generated using `gninagrid`:
```
1 CS12.48.19.binmap.gz CS12_0.48.18.binmap.gz
0 CS12.48.19.binmap.gz CS12_1.48.18.binmap.gz
0 CS12.48.19.binmap.gz CS12_2.48.18.binmap.gz
0 CS12.48.19.binmap.gz CS12_3.48.18.binmap.gz
```

As an example, imagine we want to incorporate three additional grids, `cdk_gist-dipole-dens.dx`, `cdk_gist-dipolex-dens.dx`, and `cdk_gist-gO.dx` into the input.
We would run `gninagrid`:
```
gninagrid  -r rec.pdb -l CDK2_CS12_docked.sdf.gz -g cdk_gist-dipole-dens.dx -g cdk_gist-dipolex-dens.dx -g cdk_gist-gO.dx -o CS12 --separate
```
Since `--separate` is passed, this will produce separate receptor (which includes the user provided grids) and ligand files:
```
-rw-rw-r-- 1 dkoes dkoes 8404992 Apr 28 12:55 CS12.48.19.binmap
-rw-rw-r-- 1 dkoes dkoes 7962624 Apr 28 12:55 CS12_0.48.18.binmap
-rw-rw-r-- 1 dkoes dkoes 7962624 Apr 28 12:55 CS12_1.48.18.binmap
...
```
The receptor file has 16 channels for the regular protein atom types and 3 for the provided grids.  The grid dimensions, resolution, and positioning is determined from the provided grids (which must all match).  To save (a lot of) space, the binmap files can be gzipped:
```
gzip *.binmap
```

Note that it is up to the user to ensure that the dimensions (including _total_ number of channels) of the input files match the specified dimensions in NGridLayer.

