nextflow.enable.dsl=2

process DOWNLOAD_SRA {

    tag "$srr_id"
    container 'pierrejeangouze/sra_toolkit'
    errorStrategy 'retry'
    maxRetries 2

    input:
    val srr_id

    output:
    path "${srr_id}_1.fastq.gz"
    path "${srr_id}_2.fastq.gz", optional: true

    script:
    """
    fasterq-dump --threads 12 --split-files --progress ${srr_id}
    gzip *.fastq
    """
}
    

workflow {

    // Liste definie directement ici
    srr_list = [
        'SRR10379723'
    ]

    // Afficher la liste utilisee
    println "Liste des echantillons : ${srr_list.join(', ')}"

    // Creer un canal qui emet chaque SRR separement
    srr_ch = Channel.from(srr_list)

    // Lancer une tache pour chaque SRR
    DOWNLOAD_SRA(srr_ch)
}
