FROM ubuntu:20.04

LABEL maintainer="pierrejeangouze"
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        wget \
        build-essential \
        tar \
        zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV SUBREAD_VERSION=1.4.6-p3

RUN wget https://downloads.sourceforge.net/project/subread/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-source.tar.gz && \
    tar -xzf subread-${SUBREAD_VERSION}-source.tar.gz && \
    cd subread-${SUBREAD_VERSION}-source/src && \
    make -f Makefile.Linux && \
    mv ../bin/* /usr/local/bin/ && \
    cd / && rm -rf subread-${SUBREAD_VERSION}-source*

ENTRYPOINT ["featureCounts"]

