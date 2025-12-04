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
## Prérequis : 

- Docker
- Git
- Architecture Linux/AMD64 pour les CPU
- Nextflow

### Installation de Nextflow 
Commandes bash à rentrer si Nextflow n'est pas présent sur la machine : 
```bash
sudo apt update
sudo apt install -y openjdk-21-jdk
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
nextflow -version
```

### Fichiers nécessaires
- **`main.nf`**, **`nextflow.config`**, **`Data/`**

### Lancer le pipeline

```bash
nextflow run main.nf
```

### Configuration testée : 
VM sur le Cloub IFB : 
- 8 vCPU
- 32 Go de RAM
- 200 Go de stockage

Temps d'execution : 

