FROM python:3.8

WORKDIR /rasputin

# System dependencies
RUN apt update && apt install -y \
  build-essential \
  cmake \
  libgmp-dev \
  libmpfr-dev \
  libblas3 \
  libblas-dev \
  liblapack3 \
  liblapack-dev \
  libatlas3-base \
  libatlas-base-dev \
  libgeos-dev \
  python3-tk && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## cmake
RUN wget https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.tar.gz \
    && tar xvf cmake-3.12.4-Linux-x86_64.tar.gz \
    && rm cmake-3.12.4-Linux-x86_64.tar.gz
ENV PATH="/rasputin/cmake-3.12.4-Linux-x86_64/bin:${PATH}"


## boost
RUN wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.bz2 && \
  tar xvf boost_1_69_0.tar.bz2 && \
  rm boost_1_69_0.tar.bz2
ENV BOOST_ROOT /rasputin/boost_1_69_0

# Python dependencies
ADD Pipfile .
ADD Pipfile.lock .
RUN pip install pipenv
RUN pipenv install -d --system

# Copy source code
ADD lib lib
RUN cd lib && \
    git clone https://github.com/pybind/pybind11.git --branch=v2.4.3 && \
    git clone https://gitlab.com/conradsnicta/armadillo-code.git --branch=9.800.x armadillo && \
    git clone https://github.com/boostorg/geometry.git --branch=boost-1.71.0 && \
    git clone https://github.com/HowardHinnant/date.git --branch=v2.4.1 && \
    git clone https://github.com/catchorg/Catch2.git --branch=v2.10.2 catch2 && \
    git clone https://github.com/CGAL/cgal.git --branch=releases/CGAL-5.0 cgal && \
    cd $HOME

ADD src src
ADD tests tests
ADD web web
ADD CMakeLists.txt .
ADD setup.py .

# Install package
RUN python setup.py install
