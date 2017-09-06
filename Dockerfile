FROM ubuntu:latest
ARG HYPHYVER=2.3.2
RUN apt-get update -q && \
    apt-get install curl cmake g++ ocl-icd-opencl-dev mpich python3 python3-scipy python3-tabulate r-base-dev -qqy && \
    curl -L https://github.com/veg/hyphy/archive/${HYPHYVER}.tar.gz -o /tmp/hyphy-${HYPHYVER}.tar.gz && \ 
    cd /tmp && tar -xf hyphy-${HYPHYVER}.tar.gz && \
    cd hyphy-${HYPHYVER} && cmake . && \
    make HYPHYMP install && \
    apt-get remove cmake g++ ocl-icd-opencl-dev mpich -qqy
WORKDIR /app
VOLUME /app
