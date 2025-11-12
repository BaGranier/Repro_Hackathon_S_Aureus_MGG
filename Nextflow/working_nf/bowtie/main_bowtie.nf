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
    ref = DOWNLOAD_REFERENCE()
    INDEX_CREATION(ref[0]) // le premier élément de sortie = reference.fasta
}
