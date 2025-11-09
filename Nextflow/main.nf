nextflow.enable.dsl=2

process DOWNLOAD_SRA {

    container 'pierrejeangouze/sra_toolkit'

    input:
    val srr_id

    output:
    path "*.fastq.gz"

    script:
    """
    echo "Telechargement de l'echantillon : ${srr_id}"
    fasterq-dump --threads 4 --progress ${srr_id}
    gzip *.fastq
    """
}

workflow {
    main:
        DOWNLOAD_SRA("SRR10379723")
}
