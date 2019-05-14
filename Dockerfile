FROM python:3.6

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
  libatlas-dev \
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
ADD src src
ADD tests tests
ADD web web
ADD CMakeLists.txt .
ADD setup.py .

# Install package
RUN python setup.py install

# Run tests
RUN pytest tests/
