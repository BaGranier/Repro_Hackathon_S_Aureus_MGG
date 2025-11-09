#!/bin/bash
set -e
echo "=== Mise à jour du système ==="
sudo apt update
echo "=== Installation de Java 21 ==="
sudo apt install -y openjdk-21-jdk
echo "=== Téléchargement de Nextflow ==="
curl -s https://get.nextflow.io | bash
echo "=== Installation de Nextflow dans /usr/local/bin ==="
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
echo "=== Vérification de Nextflow ==="
nextflow -version
echo "=== Exécution du pipeline ==="
cd "$(dirname "$0")"
nextflow run main.nf -with-docker
echo "✅ Exécution terminée avec succès."
