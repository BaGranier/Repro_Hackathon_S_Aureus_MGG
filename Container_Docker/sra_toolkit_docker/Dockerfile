FROM ubuntu:22.04

LABEL maintainer="pj.gouze@gmail.com"
LABEL description="Container Docker pour le SRA Toolkit"

# Installer les dépendances
RUN apt-get update && apt-get install -y \
    wget \
    tar \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Télécharger et installer SRA Toolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    mkdir -p /opt/sratoolkit && \
    mv sratoolkit.*/* /opt/sratoolkit/ && \
    ln -s /opt/sratoolkit/bin/* /usr/local/bin/ && \
    rm -rf sratoolkit.* sratoolkit.current-ubuntu64.tar.gz

# Définir le répertoire de travail
WORKDIR /data

# Monter un volume pour partager des fichiers
VOLUME ["/data"]

# Commande par défaut pour tester l’installation
CMD ["fastq-dump", "--version"]

