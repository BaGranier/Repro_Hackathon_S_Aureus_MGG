Ce dossier contient tous les containers docker

Chaque dossier contient son dockerfile permettant la construction de l'image, ainsi que les packages associés permettant d'introduire toutes les bibliothèques nécessaires à la construction de l'image. 


| **Dossier / outil** | **Description**                                                                                                        |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| **`Bowtie`**        | Contient les fichiers liés à l’alignement des lectures avec **Bowtie2** (scripts, index, paramètres).                           |
| **`FeatureCounts`** | Contient les fichiers ou scripts liés à la quantification des lectures avec **featureCounts**.                                  |
| **`Cutadapt`**      | Contient les fichiers liés au trimming/adaptation des séquences avec **Cutadapt** (paramètres, scripts).                        |
| **`TrimGalore`**    | Contient les fichiers relatifs à **TrimGalore**, wrapper autour de Cutadapt (options, logs, éventuellement scripts).            |
| **`R_Deseq2`**      | Contient le script R et les ressources nécessaires à l’analyse **DESeq2** (normalisation, tests différentiels, MA-plots, etc.). |
| **`SRA_Toolkit`**   | Contient les éléments nécessaires à l’utilisation de **SRA Toolkit** (scripts de téléchargement, métadonnées, logs).            |


