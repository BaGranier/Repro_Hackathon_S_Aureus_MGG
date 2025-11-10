FROM --platform=linux/amd64 ubuntu:14.04

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /opt

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential \
      wget \
      unzip \
      zlib1g-dev \
      ca-certificates \
      gcc-4.4 \
      g++-4.4 \
      && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget -O bowtie-0.12.7-src.zip \
      "https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download" && \
    unzip bowtie-0.12.7-src.zip && \
    rm bowtie-0.12.7-src.zip && \
    cd bowtie-0.12.7 && \
    make CC=gcc-4.4 CXX=g++-4.4 -j$(nproc) && \
    ln -sf /opt/bowtie-0.12.7/bowtie /usr/local/bin/bowtie

WORKDIR /data
ENTRYPOINT ["bowtie"]

