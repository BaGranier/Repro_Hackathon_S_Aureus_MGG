# Utiliser Ubuntu 20.04 comme base
FROM ubuntu:20.04

# Éviter les interactions lors des installations
ENV DEBIAN_FRONTEND=noninteractive

# Mettre à jour et installer les dépendances de base
RUN apt-get update && apt-get install -y \
    wget \
    python2.7 \
    python2.7-dev \
    python3 \
    python3-pip \
    build-essential \
    git \
    curl \
    ca-certificates \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Installer pip pour Python 2
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py && \
    python2.7 get-pip.py && \
    rm get-pip.py

# Créer un environnement virtuel Python 2 pour Cutadapt 1.11
RUN python2.7 -m pip install --upgrade pip setuptools cython

# Installer wget (déjà fait dans apt-get mais on peut forcer si nécessaire)
RUN apt-get update && apt-get install -y wget

# Installer Cutadapt 1.11
WORKDIR /opt
RUN wget https://github.com/marcelm/cutadapt/archive/refs/tags/v1.11.tar.gz && \
    tar -xzf v1.11.tar.gz && \
    cd cutadapt-1.11 && \
    python2.7 setup.py install

# Installer Trim Galore 0.6.6
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.6.tar.gz -O trim_galore-0.6.6.tar.gz && \
    tar -xzf trim_galore-0.6.6.tar.gz && \
    cd TrimGalore-0.6.6 && \
    chmod +x trim_galore && \
    ln -s $(pwd)/trim_galore /usr/local/bin/trim_galore

# Définir le répertoire de travail
WORKDIR /data

# Commande par défaut (exemple, peut être remplacée par ton script)
CMD ["bash"]

