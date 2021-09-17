FROM ubuntu:20.04

WORKDIR /root

# Fetch packages Dependancies
RUN apt-get update &&\
    apt-get upgrade -y &&\
    apt-get install -y\
    apt-utils 

RUN DEBIAN_FRONTEND="noninteractive" TZ="Europe/Amsterdam"  apt-get -y install tzdata


RUN apt-get install -y\
    build-essential \
    wget \
    libgmp-dev \
    libmpfr-dev \
    libblas3 \
    libboost-all-dev \
    liblapack3 \
    liblapack-dev \
    libatlas3-base \
    libatlas-base-dev \
    python3.6 \
    git \
    pybind11-dev \
    libarmadillo-dev \
    libcgal-dev \
    python3-pip \
    libcurl4 \
    libssl-dev \
    libcurl4-openssl-dev &&\
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Python Depndancies
RUN pip3 install numpy pyproj Pillow h5py lxml shapely descartes meshio pytest setuptools pipenv pyyaml requests ipython xarray

# Cmake
RUN wget https://github.com/Kitware/CMake/releases/download/v3.16.5/cmake-3.16.5.tar.gz \
    && tar -zxvf cmake-3.16.5.tar.gz \
    && rm cmake-3.16.5.tar.gz \
    && cd cmake-3.16.5 \
    && ./bootstrap \
    && make \
    && make install

# Github resources
RUN mkdir lib && cd lib && \
    git clone https://github.com/boostorg/geometry.git --branch=boost-1.71.0 &&\
    git clone https://github.com/HowardHinnant/date.git &&\
    git clone https://github.com/catchorg/Catch2.git --branch=v2.10.2 catch2 &&\
    git clone https://github.com/pybind/pybind11.git --branch=v2.4.3 &&\
    git clone https://gitlab.com/conradsnicta/armadillo-code.git --branch=9.800.x armadillo &&\
    git clone https://github.com/CGAL/cgal.git --branch=releases/CGAL-5.0 cgal

# environment variables
ENV CXX=/usr/bin/g++-9 

ENV CMAKE_MODULE_PATH=/root/lib

ENV Catch2_DIR=/root/lib/Catch2/CMake/build 

ENV date_DIR=/root/lib/date/cmake

ENV RASPUTIN_DATA_DIR=/root/rasputin_workspace/rasputin_data

# rasputin
#RUN mkdir rasputin && cd rasputin &&\ 
#    git clone https://github.com/expertanalytics/rasputin.git

#RUN cd rasputin &&\
#    pip3 install .

RUN echo "image is prepared to  install rasputin"
