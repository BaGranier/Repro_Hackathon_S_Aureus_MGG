nextflow.enable.dsl=2

process DOWNLOAD_SRA {

    tag "$srr_id"
    container 'pierrejeangouze/sra_toolkit'

    input:
    val srr_id

    output:
    path "*.fastq.gz"

    script:
    """
    echo "Telechargement de l'echantillon : ${srr_id}"
    echo "----------------------------------------------"
    fasterq-dump --threads 12 --progress ${srr_id}
    gzip *.fastq
    mkdir -p results/fastq
    mv *.fastq.gz results/fastq/
    echo "Telechargement termine : ${srr_id}"
    """
}

workflow {

    // Liste definie directement ici
    srr_list = [
        'SRR10379721',
        'SRR10379722',
        'SRR10379723',
        'SRR10379724',
        'SRR10379725',
        'SRR10379726'
    ]

    // Afficher la liste utilisee
    println "Liste des echantillons : ${srr_list.join(', ')}"

    // Creer un canal qui emet chaque SRR separement
    srr_ch = Channel.from(srr_list)

    // Lancer une tache pour chaque SRR
    DOWNLOAD_SRA(srr_ch)
}
