nextflow.enable.dsl=2

process GET_SRA {
    tag "$sra_id"
    publishDir params.outdir, mode: 'copy'

    input:
        val sra_id

    output:
        path "*.fastq.gz", emit: sra_fastq_files

    script:
        """
        echo "=== Téléchargement du fichier SRA : ${sra_id} ==="
        fasterq-dump --threads ${task.cpus} --progress ${sra_id}
        gzip *.fastq
        """
}



