#Bon Dockerfile

FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python2.7 python2.7-dev build-essential wget zlib1g-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://bootstrap.pypa.io/pip/2.7/get-pip.py -O /tmp/get-pip.py && \
    python2.7 /tmp/get-pip.py && \
    rm /tmp/get-pip.py

RUN python2.7 -m pip install --upgrade setuptools cython

WORKDIR /opt

RUN wget https://github.com/marcelm/cutadapt/archive/refs/tags/v1.11.tar.gz && \
    tar -xzf v1.11.tar.gz && \
    cd cutadapt-1.11 && \
    python2.7 setup.py install

WORKDIR /data

RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.6.tar.gz -O trim_galore-0.6.6.tar.gz && \
    tar -xzf trim_galore-0.6.6.tar.gz && \
    cd TrimGalore-0.6.6 && \
    chmod +x trim_galore && \
    ln -s $(pwd)/trim_galore /usr/local/bin/trim_galore

ENTRYPOINT ["/bin/bash", "-c"]
