nextflow.enable.dsl=2

process deseq2 {
    tag "Running"
    container 'vmichelet/r341_desq2'

    input:
    path counts
    file deseq_script

    output:
    path "*.png"

    script:
    """
    Rscript ${deseq_script} ${counts}
    """
}

workflow {
    deseq2(file("counts.txt"), file("data/deseq_fonc2.R") )
}
