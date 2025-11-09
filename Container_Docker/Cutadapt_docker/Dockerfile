# Étape 1 : utiliser une version ancienne de Python compatible avec Cutadapt 1.11
FROM python:3.6-slim

# Étape 2 : installer les dépendances système nécessaires
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Étape 3 : installer une version précise de Cutadapt
RUN pip install cutadapt==1.11

# Étape 4 : définir le point d’entrée
ENTRYPOINT ["cutadapt"]

