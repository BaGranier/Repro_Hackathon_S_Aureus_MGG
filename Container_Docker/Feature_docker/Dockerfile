# Utiliser une image de base Ubuntu
FROM ubuntu:20.04

# Définir le mainteneur
LABEL maintainer="pierrejeangouze"

# Éviter les invites interactives (ex: fuseaux horaires)
ENV DEBIAN_FRONTEND=noninteractive

# Mettre à jour et installer les dépendances nécessaires
RUN apt-get update && \
    apt-get install -y \
        wget \
        tar \
        build-essential \
        libcurl4-openssl-dev \
        libxml2-dev \
        libssl-dev \
        zlib1g-dev \
        python3 \
        python3-pip \
        r-base && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Télécharger et installer Subread 1.4.6-p3 (binaire précompilé)
RUN wget https://downloads.sourceforge.net/project/subread/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    tar -xzf subread-1.4.6-p3-Linux-x86_64.tar.gz && \
    # ✅ Correction ici : ajout de l’option -r pour copier le dossier utilities
    cp -r subread-1.4.6-p3-Linux-x86_64/bin/* /usr/local/bin/ && \
    rm -rf subread-1.4.6-p3-Linux-x86_64*

# Installer R et les packages nécessaires pour DESeq2
RUN R -e "install.packages('BiocManager', repos='http://cran.r-project.org')" && \
    R -e "BiocManager::install('DESeq2', ask=FALSE)" && \
    R -e "install.packages('optparse', repos='http://cran.r-project.org')"

# Définir le répertoire de travail
WORKDIR /data

# Commande par défaut
CMD ["bash"]

