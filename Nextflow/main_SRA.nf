process DOWNLOAD_SRA {

    // Utilisation du conteneur local
    container 'local/sra_toolkit_docker'

    input:
        val accession_id

    output:
        path "${accession_id}_1.fastq"
        path "${accession_id}_2.fastq"

    script:
    """
    echo "Téléchargement de $accession_id depuis NCBI..." 
    fasterq-dump $accession_id --split-files -O .
    """
}

workflow {
    Channel
        .from(['SRR10379726'])
        .set { sra_ids }

    DOWNLOAD_SRA(sra_ids)
}
