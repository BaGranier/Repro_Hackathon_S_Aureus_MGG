nextflow.enable.dsl=2

process TRIM_ALL_READS {
    tag "$fastq_file"
    container "bagranier/trim_galore_cutadapt"

    input:
    path fastq_file

    output:
    path "*", optional: true

    script:
    """
    echo "Processing file: ${fastq_file}"
    trim_galore -q 20 --phred33 --length 25 ${fastq_file}
    echo "Trimming complete for: ${fastq_file}"
    """
}

workflow {
    fastq_ch = Channel.fromPath('test*.fastq.gz')
    fastq_ch.view { "Detected FASTQ file: ${it}" }
    TRIM_ALL_READS(fastq_ch)
}
