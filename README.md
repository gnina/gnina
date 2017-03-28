gnina (pronounced NEE-na) is a fork of smina, which is a fork of AutoDock Vina.

gnina is **not** recommended for production use (*yet*) in molecular modeling tasks.  However, it *is* suitable as a platform for researching structure-based deep learning approaches as described in [our paper](https://arxiv.org/abs/1612.02751).  

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

If you are building for systems with different GPUs, include `-DCUDA_ARCH_NAME=All`.  
