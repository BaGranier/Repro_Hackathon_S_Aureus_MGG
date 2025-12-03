FROM ubuntu:22.04

LABEL maintainer="vmichelet"
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        wget \
        build-essential \
        tar \
        zlib1g-dev \
        ca-certificates && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENV SUBREAD_VERSION=2.0.6

WORKDIR /opt

RUN wget https://downloads.sourceforge.net/project/subread/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-source.tar.gz && \
    tar -xzf subread-${SUBREAD_VERSION}-source.tar.gz && \
    cd subread-${SUBREAD_VERSION}-source/src && \
    make -f Makefile.Linux && \
    mv ../bin/* /usr/local/bin/ && \
    cd /opt && rm -rf subread-${SUBREAD_VERSION}-source subread-${SUBREAD_VERSION}-source.tar.gz

WORKDIR /data

ENTRYPOINT ["featureCounts"]

