FROM ubuntu:22.04

WORKDIR /root

ENV DEBIAN_FRONTEND=noninteractive

# Fetch packages Dependencies
RUN apt-get update &&\
    apt-get upgrade -y &&\
    apt-get install -y\
    build-essential \
    wget \
    apt-utils \
    libgmp-dev \
    libmpfr-dev \
    libblas3 \
    libboost-all-dev \
    liblapack3 \
    liblapack-dev \
    libatlas3-base \
    libatlas-base-dev \
    python3-pip \
    git \
    pybind11-dev \
    libarmadillo-dev \
    libcgal-dev \
    python3-pip \
    libcurl4 \
    libssl-dev \
    libcurl4-openssl-dev \
    libtinfo5 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Python Dependencies
RUN python3 -m pip install numpy pyproj Pillow h5py lxml shapely descartes meshio pytest setuptools pipenv pyyaml requests

# Cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.25.2/cmake-3.25.2-linux-x86_64.tar.gz \
    && tar -zxf cmake-3.25.2-linux-x86_64.tar.gz \
    && rm cmake-3.25.2-linux-x86_64.tar.gz \
    && cp -R cmake-3.25.2-linux-x86_64/* /usr/ \
    && rm -r cmake-3.25.2-linux-x86_64

# clang+llvm 16
RUN wget https://github.com/llvm/llvm-project/releases/download/llvmorg-16.0.0/clang+llvm-16.0.0-x86_64-linux-gnu-ubuntu-18.04.tar.xz \
    && tar -xf clang+llvm-16.0.0-x86_64-linux-gnu-ubuntu-18.04.tar.xz \
    && rm clang+llvm-16.0.0-x86_64-linux-gnu-ubuntu-18.04.tar.xz \
    && cp -R clang+llvm-16.0.0-x86_64-linux-gnu-ubuntu-18.04/* /usr/ \
    && rm -r clang+llvm-16.0.0-x86_64-linux-gnu-ubuntu-18.04

ENV CXX=/usr/bin/clang++

# Github resources
RUN mkdir lib && cd lib && \
    git clone https://github.com/boostorg/geometry.git --branch=boost-1.81.0 --depth=1 &&\
    git clone https://github.com/catchorg/Catch2.git --branch=v3.3.1 --depth=1
