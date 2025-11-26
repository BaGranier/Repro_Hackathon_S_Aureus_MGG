nextflow.enable.dsl=2

process DOWNLOAD_SRA {

    tag "$srr_id"
    container 'pierrejeangouze/sra_toolkit'
    errorStrategy 'retry'
    maxRetries 2

    input:
    val srr_id

    output:
    path "*.fastq.gz"

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

    script: 
    """
    trim_galore -q 20 --phred33 --length 25 ${fastq_file}
    """
}

//--------------------------------------------------------------
// 1. DOWNLOAD REFERENCE
//--------------------------------------------------------------
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



//--------------------------------------------------------------
// 2. BOWTIE INDEX CREATION (all index files in a folder)
//--------------------------------------------------------------
process INDEX_CREATION {
    tag "index"
    container 'vmichelet/bowtie:0.12.7-samtools'
    containerOptions '--entrypoint ""'

    input:
    path fasta_file

    output:
    path "ref_index", emit: index_dir

    script:
    """
    mkdir ref_index
    bowtie-build ${fasta_file} ref_index/ref_index
    """
}



//--------------------------------------------------------------
// 3. ALIGN READS USING PROCESS SUBSTITUTION + PIPED SORT
//--------------------------------------------------------------
process ALIGN {
    tag "align"
    container 'vmichelet/bowtie:0.12.7-samtools'
    containerOptions '--entrypoint ""'
    cpus 4

    input:
    path fastq_file
    path index_dir

    output:
    path "*.bam"

    script:
    """
    bowtie -p ${task.cpus} -S ${index_dir}/ref_index <(gunzip -c ${fastq_file}) \
        | samtools sort -@ ${task.cpus} -o ${fastq_file.simpleName}.bam

    samtools index ${fastq_file.simpleName}.bam
    """
}

process FEATURECOUNTS {

    container 'pierrejeangouze/featurecounts:1.4.6-p3'

    input:
        path bam_list
        path gff

    output:
        path "counts.txt"

    script:
    """
    featureCounts \
        -t gene \
        -g ID \
        -F GFF \
        -T 4 \
        -a $gff \
        -o counts.txt ${bam_list.join(' ')}
    """
}

process deseq2 {
    tag "Running"
    container 'vmichelet/r341_desq2'

    input:
    path counts
    path reference
    path tsv
    file deseq_script

    output:
    path "*_simple.pdf"
    path "*_annotated.pdf"
    path "DESeq2_results_clean.csv"

    script:
    """
    Rscript ${deseq_script} ${counts} ${reference} ${tsv}
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
    fastq_files = DOWNLOAD_SRA(srr_ch)
    (fasta, gff) = DOWNLOAD_REFERENCE()
    index = INDEX_CREATION(fasta)
    fastq_ch = TRIM_ALL_READS(fastq_files)
    bam = ALIGN(fastq_ch, index)
    counts = FEATURECOUNTS(bam.collect(), gff)
    deseq2(counts, gff, file("data/GeneSpecificInformation_COL.tsv"), file("data/deseq.R") )
}
