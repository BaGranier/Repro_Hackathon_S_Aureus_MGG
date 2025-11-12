# ===========================
#  Bowtie 0.12.7 - Dockerfile
#  Correction : inclut bowtie-build et bowtie-inspect,
#  supprime l'ENTRYPOINT forcé
# ===========================

FROM --platform=linux/amd64 ubuntu:14.04

LABEL maintainer="bagranier"
LABEL description="Bowtie 0.12.7 with bowtie-build and bowtie-inspect enabled"
LABEL version="0.12.7"

# Empêcher les prompts interactifs
ENV DEBIAN_FRONTEND=noninteractive

# Définir le répertoire de travail
WORKDIR /opt

# Installer les dépendances minimales
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        build-essential \
        wget \
        unzip \
        zlib1g-dev \
        ca-certificates \
        gcc-4.4 \
        g++-4.4 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Télécharger et compiler Bowtie 0.12.7 depuis les sources
RUN wget -O bowtie-0.12.7-src.zip \
        "https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download" \
    && unzip bowtie-0.12.7-src.zip \
    && rm bowtie-0.12.7-src.zip \
    && cd bowtie-0.12.7 \
    && make CC=gcc-4.4 CXX=g++-4.4 -j$(nproc) \
    && ln -sf /opt/bowtie-0.12.7/bowtie /usr/local/bin/bowtie \
    && ln -sf /opt/bowtie-0.12.7/bowtie-build /usr/local/bin/bowtie-build \
    && ln -sf /opt/bowtie-0.12.7/bowtie-inspect /usr/local/bin/bowtie-inspect

# Répertoire de travail par défaut pour les utilisateurs
WORKDIR /data

# Retirer l’ENTRYPOINT problématique (on veut pouvoir lancer bowtie-build etc.)
ENTRYPOINT ["/bin/bash"]
