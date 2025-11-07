#!/bin/bash

# Nombre de cœurs à utiliser
CPUS=4

# Index Bowtie
INDEX="ref_index"

# Dossier contenant les fichiers trimmed
TRIMMED_DIR="Fichiers trimmed"

# Boucle sur tous les fichiers *_trimmed.fq.gz
for fq in "$TRIMMED_DIR"/SRR*_trimmed.fq.gz
do
    # Récupérer le nom de base sans chemin ni extension
    base=$(basename "$fq" _trimmed.fq.gz)

    echo "Traitement de $base ..."

    # Alignement + tri
    ./bowtie-0.12.7/bowtie -p $CPUS -S $INDEX <(gunzip -c "$fq") | \
    samtools sort -@ $CPUS > "${base}.bam"

    # Création de l'index BAM
    samtools index "${base}.bam"

    echo "$base terminé."
done


