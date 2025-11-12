#!/bin/bash
# Script pour forcer le démarrage de Docker, construire l'image, et exécuter Nextflow
set -e  # Arrête le script si une commande échoue
echo "=== Démarrage du service Docker ==="
# Force le démarrage selon le système
sudo systemctl start docker || sudo service docker start
echo "✅ Docker démarré."
echo "=== Construction de l'image Docker ==="
docker build -t local/sra_toolkit_docker -f Dockerfile_SRA .
echo "=== Lancement du pipeline Nextflow ==="
nextflow run main_SRA.nf -with-docker
echo "✅ Exécution terminée avec succès."
