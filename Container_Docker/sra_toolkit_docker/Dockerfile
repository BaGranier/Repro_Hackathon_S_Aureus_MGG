# Utiliser Ubuntu 22.04 comme base
FROM ubuntu:22.04

LABEL maintainer="pj.gouze@gmail.com"
LABEL description="Container Docker pour le SRA Toolkit avec Java"

# Éviter les prompts interactifs d’apt
ENV DEBIAN_FRONTEND=noninteractive

# Installer les dépendances et Java
RUN apt-get update && apt-get install -y \
    wget \
    tar \
    gzip \
    openjdk-11-jre \
    && rm -rf /var/lib/apt/lists/*

# Télécharger et installer le SRA Toolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    mkdir -p /opt/sratoolkit && \
    mv sratoolkit.*/* /opt/sratoolkit/ && \
    ln -s /opt/sratoolkit/bin/* /usr/local/bin/ && \
    rm -rf sratoolkit.* sratoolkit.current-ubuntu64.tar.gz

# Définir le répertoire de travail
WORKDIR /data

# Monter un volume pour le partage de données
VOLUME ["/data"]

# Vérifier les versions installées au démarrage
CMD ["bash", "-c", "java -version && fastq-dump --version"]

