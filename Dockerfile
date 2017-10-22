FROM ubuntu:latest
ENV LC_ALL=C.UTF-8
ARG HYPHYVER=2.3.2
RUN apt-get update -q && \
    apt-get install curl cmake g++ ocl-icd-opencl-dev mpich python3 python3-requests python3-scipy python3-tabulate python3-xlsxwriter r-base r-base-dev -qqy && \
    curl -L https://github.com/veg/hyphy/archive/${HYPHYVER}.tar.gz -o /tmp/hyphy-${HYPHYVER}.tar.gz && \ 
    cd /tmp && tar -xf hyphy-${HYPHYVER}.tar.gz && \
    cd hyphy-${HYPHYVER} && cmake . && \
    make HYPHYMP HYPHYMPI install && \
    echo "install.packages('ggplot2', repos='http://cran.us.r-project.org')" | R --no-save && \
    rm -rf /tmp/hyphy-${HYPHYVER}* && \
    apt-get remove cmake g++ ocl-icd-opencl-dev -qqy
WORKDIR /app
VOLUME /app
