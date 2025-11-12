nextflow.enable.dsl=2

process TRIM_READS {

    tag "$fastq_file"
    container 'pierrejeangouze/trimm_galore_cutadapt:1.11'

    input:
    path fastq_file

    output:
    path "*.fq.gz"

    script:
    """
    echo "Processing file: ${fastq_file}"
    trim_galore -q 20 --phred33 --length 25 ${fastq_file}
    echo "Trimming complete for: ${fastq_file}"
    """
}

workflow {
    fastq_file = file('SRR10379723.fastq.gz')
    println "Trimming file: ${fastq_file}"
    TRIM_READS(fastq_file)
}
