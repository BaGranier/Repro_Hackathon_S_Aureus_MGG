# 🔹 Changement 1 : utiliser AMD64 car r-base:3.4.1 n'existe pas pour ARM
FROM --platform=linux/amd64 r-base:3.4.1

# Mainteneur
LABEL maintainer="Pierre-Jean <pj.gouze@gmail.com>"

# 🔹 Changement 2 : réécrire sources.list pour n'utiliser que Stretch archive
#   - On ignore les clés expirées
RUN echo "deb [trusted=yes] http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list && \
    echo "deb-src [trusted=yes] http://archive.debian.org/debian stretch main contrib non-free" >> /etc/apt/sources.list && \
    echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99no-check-valid && \
    echo 'Acquire::AllowInsecureRepositories "true";' >> /etc/apt/apt.conf.d/99no-check-valid

# 🔹 Changement 3 : installation minimale des dépendances système
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev \
        libxml2-dev \
        libssl-dev \
        wget \
        ca-certificates \
        gnupg \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# 🔹 Changement 4 : installer DESeq2 via BiocManager (évite les wget manuels et les 404)
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='http://cran.r-project.org'); \
          BiocManager::install(version='3.5'); \
          BiocManager::install(c('DESeq2'))"

# 🔹 Vérification de l'installation
RUN R -e "library(DESeq2); packageVersion('DESeq2')"

# Dossier de travail
WORKDIR /data

