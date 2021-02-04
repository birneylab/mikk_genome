FROM gcc:5.4.0

MAINTAINER Ian Brettell <brettell@ebi.ac.uk>

# RUN apt-get update && apt-get -y upgrade
# RUN apt-get -y install \
#       build-essential \
#       wget \
#       bzip2 \
#       ca-certificates \
#       libglib2.0-0 \
#       libxext6 \
#       libsm6 \
#       libxrender1 \
#       git \
#       locales
# RUN rm -rf /var/lib/apt/lists/*
# RUN apt-get clean
# 
# # Set locale
#
# RUN locale-gen en_US.UTF-8

# Install `cmake`

RUN wget https://github.com/Kitware/CMake/releases/download/v3.19.4/cmake-3.19.4-Linux-x86_64.sh
RUN chmod +x ./cmake-3.19.4-Linux-x86_64.sh
RUN ./cmake-3.19.4-Linux-x86_64.sh

# Download SViper

RUN git clone --recursive https://github.com/smehringer/SViper.git
RUN cd SViper
RUN mkdir build
RUN cd build
RUN cmake ..
RUN make sviper
