# Bowtie 0.12.7 – Modern build + Samtools 1.13
# ============================================
FROM ubuntu:22.04

LABEL maintainer="bagranier"
LABEL description="Bowtie 0.12.7 rebuilt on Ubuntu 22.04 with GCC compatibility and Samtools"
LABEL version="0.12.7+samtools1.13"

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# Installer les dépendances nécessaires à la compilation + Samtools
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
        samtools && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Télécharger et compiler Bowtie
RUN wget -O bowtie-0.12.7-src.zip \
        "https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download" && \
    unzip bowtie-0.12.7-src.zip && \
    rm bowtie-0.12.7-src.zip && \
    cd bowtie-0.12.7 && \
    # ⚙️ Patch de compatibilité : forcer C++98 et désactiver les warnings fatals
    make CXXFLAGS="-std=c++98 -fpermissive -Wno-narrowing" -j$(nproc) && \
    ln -sf /opt/bowtie-0.12.7/bowtie /usr/local/bin/bowtie && \
    ln -sf /opt/bowtie-0.12.7/bowtie-build /usr/local/bin/bowtie-build && \
    ln -sf /opt/bowtie-0.12.7/bowtie-inspect /usr/local/bin/bowtie-inspect

# Répertoire de travail par défaut
WORKDIR /data

# Entrée interactive
ENTRYPOINT ["/bin/bash"]

