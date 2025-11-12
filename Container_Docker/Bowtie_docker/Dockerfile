# ===========================
# Bowtie 0.12.7 – Modern build
# ===========================
FROM ubuntu:22.04

LABEL maintainer="bagranier"
LABEL description="Bowtie 0.12.7 rebuilt on Ubuntu 22.04 with full toolset"
LABEL version="0.12.7"

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# Installer les dépendances nécessaires à la compilation
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        unzip \
        zlib1g-dev \
        ca-certificates && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Télécharger et compiler Bowtie
RUN wget -O bowtie-0.12.7-src.zip \
        "https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download" && \
    unzip bowtie-0.12.7-src.zip && \
    rm bowtie-0.12.7-src.zip && \
    cd bowtie-0.12.7 && \
    make -j$(nproc) && \
    ln -sf /opt/bowtie-0.12.7/bowtie /usr/local/bin/bowtie && \
    ln -sf /opt/bowtie-0.12.7/bowtie-build /usr/local/bin/bowtie-build && \
    ln -sf /opt/bowtie-0.12.7/bowtie-inspect /usr/local/bin/bowtie-inspect

# Répertoire de travail par défaut
WORKDIR /data

# Entrée interactive
ENTRYPOINT ["/bin/bash"]
