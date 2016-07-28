#!/bin/sh

# make -j8 -O BASE=/usr CUDA_INCLUDE=/opt/cuda/include OPENBABEL_INCLUDE=/usr/include/openbabel-2.0  NVCC=/opt/cuda/bin/nvcc LDFLAGS='-L/usr/lib -L. -L/usr/local/lib -L/opt/cuda/lib64'

if [ "$(hostname)" = "t410" ]
then
    make -j4 BASE=/usr CUDA_INCLUDE=/opt/cuda/include OPENBABEL_INCLUDE=/usr/include/openbabel-2.0  NVCC=/opt/cuda/bin/nvcc LDFLAGS='-L/usr/lib -L. -L/usr/local/lib -L/opt/cuda/lib64'
else
    make -j8 BASE=/usr CUDA_INCLUDE=/usr/local/cuda/include OPENBABEL_INCLUDE=/usr/include/openbabel-2.0  NVCC=/usr/local/cuda/bin/nvcc LDFLAGS='-L/usr/lib -L. -L/usr/local/lib -L/usr/local/cuda/lib64'\
         2>&1 | sed -r "s/\\(([0-9]+)\\)/:\\1/g" 1>&2;
fi

cp smina.gpu ../../../test/
