nextflow.enable.dsl=2

process DOWNLOAD_REFERENCE {
    tag "download_ref"

    output:
    path "reference.fasta"
    path "reference.gff"

    script:
    """
    wget -q -O reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
    wget -q -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"
    ls -lh reference.fasta reference.gff
    """
}

process INDEX_CREATION {
    tag "index_creation"
    container 'vmichelet/bowtie:0.12.7-samtools'
    containerOptions '--entrypoint ""'

    input:
    path fasta_file

    output:
    path "ref_index_folder"

    script:
    """
    mkdir -p ref_index_folder
    bowtie-build $fasta_file ref_index_folder/ref_index
    ls -lh ref_index_folder
    """
}

process ALIGN_READS {
    tag "align_reads"
    container 'vmichelet/bowtie:0.12.7-samtools'
    containerOptions '--entrypoint ""'

    input:
    tuple path(index_dir), path(reads)

    output:
    path "*.bam"
    path "*.bai"

    cpus 4

    script:
    """
    echo "Alignement du fichier: $reads"

    # Copier l’index dans le dossier de travail
    cp -r $index_dir/* ./

    # Base name pour sortie
    base=\$(basename $reads _trimmed.fq.gz)

    # Décompression temporaire
    gunzip -c "$reads" > temp_reads.fq

    # Alignement Bowtie 0.12.7
    bowtie -p $task.cpus ref_index temp_reads.fq > "\${base}.sam"

    # Conversion SAM → BAM et tri
    samtools view -bS "\${base}.sam" | samtools sort -@ $task.cpus -o "\${base}.bam"

    # Index BAM
    samtools index "\${base}.bam"

    # Nettoyage
    rm "\${base}.sam" temp_reads.fq

    ls -lh "\${base}.bam" "\${base}.bam.bai"
    """
}

workflow {

    // 1️⃣ Télécharger la référence
    ref = DOWNLOAD_REFERENCE()

    // 2️⃣ Créer l’index
    index_ch = INDEX_CREATION(ref[0])

    // 3️⃣ Channel de tous les FASTQ dans le dossier courant
    reads_ch = Channel.fromPath("*_trimmed.fq.gz")
    reads_ch.view { "FastQ détecté: $it" }

    // 4️⃣ Tuple (index, reads) dans le bon ordre
    reads_with_index = reads_ch.combine(index_ch)
        .map { tuple(it[1], it[0]) }
    reads_with_index.view { "Tuple (index, reads): $it" }

    // 5️⃣ Lancer l’alignement
    ALIGN_READS(reads_with_index)
}
