FROM ubuntu:latest
ADD https://github.com/veg/hyphy/archive/v2.3.3.tar.gz /tmp/hyphy-2.3.3.tar.gz
RUN apt-get update -q && \
    apt-get install cmake g++ ocl-icd-opencl-dev mpich python3 -qqy && \
    cd /tmp && tar -xf hyphy-2.3.3.tar.gz && \
    cd hyphy-2.3.3 && cmake . && \
    make HYPHYMP install && \
    apt-get remove cmake g++ ocl-icd-opencl-dev mpich -qqy
WORKDIR /app
VOLUME /app
