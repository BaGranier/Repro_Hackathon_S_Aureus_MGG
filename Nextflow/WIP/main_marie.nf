nextflow.enable.dsl=2

// Assurez-vous que le chemin d'inclusion est correct
include { GET_REF_GENOME } from "./processes/GET_REF_GENOME/main" 

// Définition du paramètre

workflow {
    // 1. Appel du processus
    ch_ref_genome = GET_REF_GENOME(params.ref_genome).ref_genome_file 

    // 2. Affichage du résultat pour confirmation (le plus simple)
    ch_ref_genome.view { "Référence génome téléchargée : $it" }
}