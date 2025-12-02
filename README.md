# Repro_Hackathon_S_Aureus_MGG
Valentine Michelet, Pierre-Jean Gouze et Baptiste Granier


Ce projet reproduit le pipeline **RNA-seq** de l’article *“Intracellular Staphylococcus aureus persisters upon antibiotic exposure”* de manière **reproductible**, depuis les fichiers FASTQ jusqu’aux analyses différentielles (et production des MA-plots DESeq2).

L’objectif est de garantir la reproductibilité complète des résultats : mêmes versions de packages, mêmes outils.

---

## Structure du dépôt

| Dossier / fichier | Description |
|-------------------|--------------|
| **`main.nf`** | Script **Nextflow** principal décrivant l’ensemble du pipeline. |
| **`nextflow.config`** | Fichier de configuration du pipeline : définit les images Docker à utiliser (appel depuis DockerHub). |
| **`Container_docker/`** | Contient un Dockerfile pour chaque outil utilisé (Bowtie, Cutadapt, FeatureCounts, DESeq2, etc.). Les images correspondantes sont disponibles sur DockerHub. |
| **`Data/`** | Données d’entrée :<br>• `config.csv` - table de description des échantillons (nom, URL FASTQ, réplicat, condition)<br>• script R (analyse DESeq2). |
| **`to_delete/`** | Dossier temporaire pour fichiers/données à valider avant suppression définitive. |

---

## Objectif


Reproduire les figures principales de l’article à partir des données publiques, dans un environnement **conteneurisé** et **traçable** via Nextflow et Docker.

---



### Fichiers nécessaires
- **`main.nf`**, **`nextflow.config`**, **`Data/`**

### Lancer le pipeline

```bash
nextflow run main.nf