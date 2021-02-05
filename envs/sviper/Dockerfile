FROM gcc:5.4.0

MAINTAINER Ian Brettell <brettell@ebi.ac.uk>

# Install `cmake`
RUN apt-get -y update

RUN wget https://github.com/Kitware/CMake/releases/download/v3.19.4/cmake-3.19.4.tar.gz
RUN tar -xzf cmake-3.19.4.tar.gz
RUN ./cmake-3.19.4/bootstrap && make && make install

# Install `SViper`

RUN git clone --recursive https://github.com/smehringer/SViper.git
RUN cd SViper
RUN mkdir build
RUN cd build
RUN cmake ..
RUN make sviper
