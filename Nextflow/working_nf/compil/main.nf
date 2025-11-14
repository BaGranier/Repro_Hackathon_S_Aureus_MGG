nextflow.enable.dsl=2

process DOWNLOAD_SRA {

    tag "$srr_id"
    container 'pierrejeangouze/sra_toolkit'
    errorStrategy 'retry'
    maxRetries 2

    input:
    val srr_id

    output:
    path "${srr_id}.fastq.gz"

    script:
    """
    # Telechargement
    fasterq-dump --threads 12 --split-files --progress ${srr_id}
    # Compression
    gzip *.fastq
    """
}


process TRIM_ALL_READS {
    tag "${fastq_file.simpleName}"
    container "bagranier/trim_galore_cutadapt"

    input:
    path fastq_file

    output:
    path "*.fq.gz"
    path "*report.txt", optional: true

    script:
    """
    trim_galore -q 20 --phred33 --length 25 ${fastq_file}
    """
}

process DOWNLOAD_REFERENCE {
    tag "download_ref"

    output:
    path "reference.fasta"
    path "reference.gff"

    script:
    """
    wget -q -O reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
    wget -q -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
    """
}

process INDEX_CREATION {
    tag "INDEX_CREATION"
    container 'bagranier/bowtie:0.12.7-fixed'
    containerOptions '--entrypoint ""'

    input:
    path fasta_file

    output:
    path "ref_index.*"

    script:
    """
    bowtie-build $fasta_file ref_index
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
    fastq_files = DOWNLOAD_SRA(srr_ch)
    ref = DOWNLOAD_REFERENCE()
    INDEX_CREATION(ref[0]) // le premier élément de sortie = reference.fasta
    TRIM_ALL_READS(fastq_files)
}
