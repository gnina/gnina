# Maximum Compatibility Build - DO NOT USE FOR PRODUCTION
# Build off of centos7 (ubi 7 didn't have static libc)

FROM nvidia/cuda:10.0-cudnn7-devel-centos7

RUN yum update -y

# get packages 
WORKDIR /root
RUN yum -y install openssl-devel bzip2-devel wget make xz git epel-release \
    bzip2 glibc-static dh-autoreconf curl-devel expat-devel gettext-devel \
    perl-devel zlib-devel hdf5-devel atlas-devel atlas-static \
    lapack-devel lapack-static   python3
RUN yum -y install boost169-devel boost169-test boost169-regex \
    boost169-python3-devel hdf5-devel  dh-autoreconf hdf5-static boost169-static zlib-static openblas-static

RUN yum -y install cuda-11-2

# build compatible gcc
# download and extract source 
ENV gcc_version "6.2.0"
RUN wget --no-verbose \
        http://ftpmirror.gnu.org/gcc/gcc-${gcc_version}/gcc-${gcc_version}.tar.bz2 && \
    tar xf gcc-${gcc_version}.tar.bz2
RUN wget --no-verbose \
        https://gmplib.org/download/gmp/gmp-6.1.0.tar.xz && \
    tar xf gmp-6.1.0.tar.xz && \
    mv gmp-6.1.0 gcc-${gcc_version}/gmp
RUN wget --no-verbose \
        ftp://ftp.gnu.org/gnu/mpc/mpc-1.0.3.tar.gz && \
    tar xf mpc-1.0.3.tar.gz && \
    mv mpc-1.0.3 gcc-${gcc_version}/mpc
RUN wget --no-verbose \
				https://ftp.gnu.org/gnu/mpfr/mpfr-3.1.4.tar.xz && \
    tar xf mpfr-3.1.4.tar.xz && \
    mv mpfr-3.1.4 gcc-${gcc_version}/mpfr

# compile and install gcc
RUN mkdir gcc-${gcc_version}_build && \
    cd gcc-${gcc_version}_build && \
    ../gcc-${gcc_version}/configure \
        --prefix=/usr \
        --disable-multilib \
        --enable-languages=c,c++ \
        --enable-libstdcxx-threads \
        --enable-libstdcxx-time \
        --enable-shared \
        --enable-__cxa_atexit \
        --disable-libunwind-exceptions \
        --disable-libada \
        --host x86_64-redhat-linux-gnu \
        --build x86_64-redhat-linux-gnu \
        --with-default-libstdcxx-abi=gcc4-compatible
RUN cd gcc-${gcc_version}_build && make -j8
RUN cd gcc-${gcc_version}_build && make install

# register new libraries with `ldconfig`
RUN echo "/usr/local/lib64" > usrLocalLib64.conf && \
    mv usrLocalLib64.conf /etc/ld.so.conf.d/ && \
    ldconfig

# clean out all the garbage
RUN rm -rf ~/${gcc_release} ~/{gcc_release}_build 

# symlink so nvcc uses the "right" gcc
RUN ln -s /usr/bin/gcc /usr/local/cuda/bin/gcc
RUN ln -s /usr/bin/g++ /usr/local/cuda/bin/g++


# set up acceptable cmake version
RUN wget https://github.com/Kitware/CMake/releases/download/v3.18.6/cmake-3.18.6-Linux-x86_64.tar.gz
RUN tar xzf cmake-3.18.6-Linux-x86_64.tar.gz
RUN cp cmake-3.18.6-Linux-x86_64/bin/* /usr/local/bin
RUN cp -r cmake-3.18.6-Linux-x86_64/share/cmake-3.18 /usr/local/share/

# set up eigen
RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
RUN tar xzf eigen-3.3.7.tar.gz
RUN cd eigen-3.3.7; mkdir build; cd build; cmake ..; make install

# build and install openbabel
# Need to be new enough to include this fix: https://github.com/openbabel/openbabel/pull/2328
# todo: change to release version once a release exists with this fixes
RUN git clone https://github.com/openbabel/openbabel.git; \
    cd openbabel; \
    git checkout openbabel-3-1-1; \
    mkdir build; cd build;  \
    cmake -DPTHREAD_LIBRARY=/lib64/libpthread.a  -DBUILD_SHARED=OFF -DPYTHON_BINDINGS=OFF  -DWITH_COORDGEN=0 -DWITH_MAEPARSER=0 -DRUN_SWIG=OFF -DCMAKE_CXX_FLAGS="-std=c++14 -fPIC" -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_C_COMPILER=/usr/bin/gcc ..;  \
    make -j8; make install

#for versioning need -C option of git
RUN git clone https://github.com/git/git.git; \
    cd git; git checkout v2.25.4; \
    make prefix=/usr all; make prefix=/usr install
    
# install molgrid 
RUN yum install -y python3-devel
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install numpy pytest pyquaternion -U
# build libmolgrid
RUN python3.6 -m pip install --upgrade pip
RUN python3.6 -m pip install numpy pytest pyquaternion auditwheel twine
RUN git clone https://github.com/gnina/libmolgrid.git; \
    cd libmolgrid; \
    mkdir build; cd build; \
    cmake -DBoost_INCLUDE_DIR=/usr/include/boost169/ -DBOOST_LIBRARYDIR=/usr/lib64/boost169/ ..; \
    make -j8; make install; \
    cmake -DBoost_INCLUDE_DIR=/usr/include/boost169/ -DBOOST_LIBRARYDIR=/usr/lib64/boost169/ -DBUILD_STATIC=ON -DBoost_USE_STATIC_LIBS=ON -DBUILD_SHARED=OFF ..; \
    make -j8 libmolgrid_static; cp lib64/libmolgrid.a /usr/local/lib/libmolgrid.a 
    
RUN git clone https://github.com/gflags/gflags.git; \
    cd gflags; git checkout v2.2.2; \
    mkdir build; cd build; \
    cmake ..; make -j8; \
    make install

RUN git clone https://github.com/google/glog.git; \
    cd glog; git checkout v0.4.0; \
    mkdir build; cd build; \
    cmake ..  -DHAVE_PWD_H=0 ; make -j8; \
    make install
    
RUN git clone https://github.com/protocolbuffers/protobuf.git; \
    cd protobuf; git checkout v3.6.1.3; \
    git submodule update --init --recursive; \
    ./autogen.sh; ./configure; \
    make -j8; \
    make install
    
RUN git clone  https://github.com/erget/libaec.git; \
    cd libaec; git checkout v0.3.4; \
    mkdir build; cd build; \
    cmake -DBUILD_SHARED_LIBS=OFF ..; make -j8; \
    make install    
                
RUN git clone https://github.com/gnina/gnina.git; \
    cd gnina; \
    git checkout v1.0.1; \
     mkdir build; cd build; \
    cmake .. -DAtlas_CBLAS_LIBRARY=/usr/lib64/libopenblas.a  -DBUILD_SHARED_LIBS=OFF -DBoost_INCLUDE_DIR=/usr/include/boost169/ \
    -DPYTHON_EXECUTABLE=/usr/bin/python3 -DBOOST_LIBRARYDIR=/usr/lib64/boost169/ -DHAVE_PYTHON=0 -DBUILD_python=OFF \
    -DLIBMOLGRID_LIBRARY=/usr/local/lib/libmolgrid.a -DCUDA_ARCH_NAME=All -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11; \
    make -j8 gnina; 
    
#make "mostly static" executable
RUN cd gnina/build/gninasrc; \
/usr/bin/c++  -Wno-deprecated-declarations -Wno-unknown-pragmas -O3 -DNDEBUG   CMakeFiles/gnina.dir/main/main.cpp.o CMakeFiles/gnina.dir/version.cpp.o CMakeFiles/gnina.dir/lib/atom_constants.cpp.o CMakeFiles/gnina.dir/lib/bfgs.cu.o CMakeFiles/gnina.dir/lib/box.cpp.o CMakeFiles/gnina.dir/lib/builtinscoring.cpp.o CMakeFiles/gnina.dir/lib/cache.cpp.o CMakeFiles/gnina.dir/lib/cache_gpu.cpp.o CMakeFiles/gnina.dir/lib/cnn_scorer.cpp.o CMakeFiles/gnina.dir/lib/cnn_data.cpp.o CMakeFiles/gnina.dir/lib/coords.cpp.o CMakeFiles/gnina.dir/lib/custom_terms.cpp.o CMakeFiles/gnina.dir/lib/device_buffer.cpp.o CMakeFiles/gnina.dir/lib/everything.cpp.o CMakeFiles/gnina.dir/lib/flexinfo.cpp.o CMakeFiles/gnina.dir/lib/GninaConverter.cpp.o CMakeFiles/gnina.dir/lib/grid.cpp.o CMakeFiles/gnina.dir/lib/grid_gpu.cu.o CMakeFiles/gnina.dir/lib/model.cpp.o CMakeFiles/gnina.dir/lib/molgetter.cpp.o CMakeFiles/gnina.dir/lib/monte_carlo.cpp.o CMakeFiles/gnina.dir/lib/mutate.cpp.o CMakeFiles/gnina.dir/lib/my_pid.cpp.o CMakeFiles/gnina.dir/lib/naive_non_cache.cpp.o CMakeFiles/gnina.dir/lib/non_cache.cpp.o CMakeFiles/gnina.dir/lib/non_cache_cnn.cpp.o CMakeFiles/gnina.dir/lib/obmolopener.cpp.o CMakeFiles/gnina.dir/lib/parallel_mc.cpp.o CMakeFiles/gnina.dir/lib/parallel_progress.cpp.o CMakeFiles/gnina.dir/lib/parse_pdbqt.cpp.o CMakeFiles/gnina.dir/lib/pdb.cpp.o CMakeFiles/gnina.dir/lib/PDBQTUtilities.cpp.o CMakeFiles/gnina.dir/lib/quasi_newton.cpp.o CMakeFiles/gnina.dir/lib/quaternion.cu.o CMakeFiles/gnina.dir/lib/random.cpp.o CMakeFiles/gnina.dir/lib/result_info.cpp.o CMakeFiles/gnina.dir/lib/ssd.cpp.o CMakeFiles/gnina.dir/lib/szv_grid.cpp.o CMakeFiles/gnina.dir/lib/terms.cpp.o CMakeFiles/gnina.dir/lib/weighted_terms.cpp.o CMakeFiles/gnina.dir/lib/conf.cpp.o CMakeFiles/gnina.dir/lib/conf_gpu.cu.o CMakeFiles/gnina.dir/lib/gpucode.cu.o CMakeFiles/gnina.dir/lib/model.cu.o CMakeFiles/gnina.dir/lib/non_cache_gpu.cu.o CMakeFiles/gnina.dir/lib/precalculate_gpu.cu.o CMakeFiles/gnina.dir/lib/tree_gpu.cu.o CMakeFiles/gnina.dir/lib/user_opts.cpp.o CMakeFiles/gnina.dir/cmake_device_link.o  -o ../bin/gnina  \
 -L/usr/local/cuda/targets/x86_64-linux/lib  -Wl,-rpath,/usr/lib64/boost169:/usr/local/lib:/usr/local/cuda/lib64:/usr/lib64/atlas: \
   -Wl,--whole-archive ../caffe/lib/libcaffe.a -Wl,--no-whole-archive  \
  /usr/lib64/boost169/libboost_program_options.a \
  /usr/lib64/boost169/libboost_system.a \
  /usr/lib64/boost169/libboost_iostreams.a \
  /usr/lib64/boost169/libboost_timer.a \
  /usr/lib64/boost169/libboost_thread.a \
  /usr/lib64/boost169/libboost_serialization.a \
  /usr/lib64/boost169/libboost_filesystem.a \
  /usr/lib64/boost169/libboost_date_time.a \
  /usr/lib64/boost169/libboost_regex.a \
  /usr/lib64/boost169/libboost_unit_test_framework.a \
  /usr/lib64/boost169/libboost_chrono.a \
  /usr/lib64/boost169/libboost_atomic.a \
  /usr/local/lib/libopenbabel.a \
  /usr/local/lib/libmolgrid.a \
  ../caffe/lib/libcaffeproto.a \
  /usr/local/lib64/libglog.a \
  /usr/local/lib/libgflags.a \
  /usr/local/lib/libprotobuf.a \
  /usr/lib64/libdl.a \
  /usr/lib64/libhdf5_hl_cpp.a \
  /usr/lib64/libhdf5_hl.a  \
  /usr/lib64/libhdf5_cpp.a \
  /usr/lib64/libhdf5.a \
  /usr/lib64/libz.a \
  /usr/local/cuda/lib64/libcurand_static.a  \
  /usr/lib64/libcudnn_static.a \
  /usr/local/cuda/lib64/libcublas_static.a \
  /usr/local/cuda/lib64/libculibos.a \
  /usr/lib64/liblapack.a \
  /usr/lib64/libopenblas64.a  \
  /usr/lib64/atlas/libatlas.a \
  /usr/local/lib/libsz.a \
  /usr/local/lib/libaec.a \
  /usr/local/cuda/lib64/libcudart_static.a  \
   /usr/local/cuda/targets/x86_64-linux/lib/libcudadevrt.a \
   /usr/local/cuda/lib64/libcublasLt_static.a \
     -lm -lpthread -ldl -lrt 

#./bin/gnina -r ../test/gnina/data/184l_rec.pdb -l ../test/gnina/data/184l_lig.sdf --autobox_ligand ../test/gnina/data/184l_lig.sdf --seed 0


