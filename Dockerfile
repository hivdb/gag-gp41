FROM ubuntu:latest
ARG HYPHYVER=2.2.7
RUN apt-get update -q && \
    apt-get install curl cmake g++ ocl-icd-opencl-dev mpich python3 -qqy && \
    curl -L https://github.com/veg/hyphy/archive/${HYPHYVER}.tar.gz -o /tmp/hyphy-${HYPHYVER}.tar.gz && \ 
    cd /tmp && tar -xf hyphy-${HYPHYVER}.tar.gz && \
    cd hyphy-${HYPHYVER} && cmake . && \
    make HYPHYMP install && \
    apt-get remove cmake g++ ocl-icd-opencl-dev mpich -qqy
WORKDIR /app
VOLUME /app
