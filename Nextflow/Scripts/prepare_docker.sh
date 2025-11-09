#!/usr/bin/env bash

set -e  # stop si erreur
echo "🚀 Initialisation du build Docker SRA Toolkit"
echo "---------------------------------------------"

# === Vérification Docker ===
if ! command -v docker &> /dev/null; then
    echo "❌ Docker n'est pas disponible. Lance Docker Desktop puis réessaie."
    exit 1
fi

echo "✅ Docker détecté : $(docker --version)"
echo "Vérification du daemon..."
if ! docker info &> /dev/null; then
    echo "❌ Docker daemon inactif — ouvre Docker Desktop sous Windows."
    exit 1
fi
echo "✅ Docker daemon actif"

# === Téléchargement image de base ===
echo "➡️  Téléchargement de l'image de base ubuntu:22.04"
docker pull ubuntu:22.04

# === Construction image locale ===
echo "➡️  Construction de l'image locale 'local/sra_toolkit_docker'"
docker build --no-cache -t local/sra_toolkit_docker -f Dockerfile_SRA .

# === Vérification Java et SRA Toolkit ===
echo "➡️  Vérification du contenu de l'image :"
docker run --rm local/sra_toolkit_docker java -version || echo "⚠️  Java non trouvé"
docker run --rm local/sra_toolkit_docker fastq-dump --version || echo "⚠️  SRA Toolkit non trouvé"

echo "---------------------------------------------"
echo "✅ Image 'local/sra_toolkit_docker' prête à l'emploi !"
echo "   → Utilisable avec Nextflow :"
echo "     nextflow run main_SRA.nf -with-docker"
echo "---------------------------------------------"
