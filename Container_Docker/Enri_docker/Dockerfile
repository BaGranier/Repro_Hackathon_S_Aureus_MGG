# -----------------------------
# Dockerfile pour EnrichmentBrowser 2.14.3
# Compatible R 3.1.3 et Bioconductor 2.14
# -----------------------------

FROM --platform=linux/amd64 ubuntu:16.04

# Mettre à jour les dépôts et installer les dépendances système
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        libcurl4-openssl-dev \
        libxml2-dev \
        libssl-dev \
        libreadline-dev \
        libx11-dev \
        libxt-dev \
        libpng-dev \
        libjpeg-dev \
        libcairo2-dev \
        libbz2-dev \
        zlib1g-dev \
        liblzma-dev \
        wget \
        unzip \
        git \
        ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Télécharger et compiler R 3.1.3
WORKDIR /tmp
RUN wget https://cran.r-project.org/src/base/R-3/R-3.1.3.tar.gz && \
    tar -xzf R-3.1.3.tar.gz && \
    cd R-3.1.3 && \
    ./configure --enable-R-shlib --with-readline --with-x && \
    make && \
    make install

# Ajouter R au PATH
ENV PATH="/usr/local/bin:${PATH}"

# Installer Bioconductor 2.14 et EnrichmentBrowser 2.14.3
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite(version='2.14')"
RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite('EnrichmentBrowser')"

# Définir le dossier de travail
WORKDIR /home/docker

# Commande par défaut
CMD ["R"]

