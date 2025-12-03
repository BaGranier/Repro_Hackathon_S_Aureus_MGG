FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-dev build-essential wget zlib1g-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Mettre à jour pip / setuptools / cython
RUN python3 -m pip install --upgrade pip setuptools cython

WORKDIR /opt

# Installer Cutadapt, version actuelle
RUN python3 -m pip install cutadapt

# Installer TrimGalore!, version actuelle
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.tar.gz -O /tmp/trim_galore-0.6.10.tar.gz && \
    tar -xzf /tmp/trim_galore-0.6.10.tar.gz -C /opt && \
    ln -s /opt/TrimGalore-0.6.10/trim_galore /usr/local/bin/trim_galore && \
    chmod +x /usr/local/bin/trim_galore

WORKDIR /data

ENTRYPOINT ["/bin/bash", "-c"]

