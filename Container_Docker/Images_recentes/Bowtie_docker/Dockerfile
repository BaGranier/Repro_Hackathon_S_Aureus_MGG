# Bowtie 1.3.1 – Modern build + Samtools 1.19
# ==========================================
FROM ubuntu:24.04

LABEL maintainer="vmichelet"
LABEL description="Bowtie 1.3.1 rebuilt on Ubuntu 24.04 with Samtools"
LABEL version="1.3.1+samtools1.19"

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# Installer dépendances pour compilation
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        unzip \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        ca-certificates \
        samtools \
        python3 \
        python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Télécharger et compiler Bowtie 1.3.1
RUN wget https://github.com/BenLangmead/bowtie/archive/refs/tags/v1.3.1.zip -O bowtie-1.3.1.zip && \
    unzip bowtie-1.3.1.zip && \
    rm bowtie-1.3.1.zip && \
    cd bowtie-1.3.1 && \
    make -j$(nproc) && \
    ln -sf /opt/bowtie-1.3.1/bowtie /usr/local/bin/bowtie && \
    ln -sf /opt/bowtie-1.3.1/bowtie-build /usr/local/bin/bowtie-build && \
    ln -sf /opt/bowtie-1.3.1/bowtie-inspect /usr/local/bin/bowtie-inspect

# Répertoire de travail par défaut
WORKDIR /data

# Entrée interactive
ENTRYPOINT ["/bin/bash"]

