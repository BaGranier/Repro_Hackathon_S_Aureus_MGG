nextflow.enable.dsl=2

process deseq2 {
    tag "Running"
    container 'vmichelet/r341_desq2'

    input:
    path counts
    path aureo
    file deseq_script

    output:
    path "*.png"

    script:
    """
    Rscript ${deseq_script} ${counts} ${aureo}
    """
}

workflow {
    deseq2(file("counts.txt"), file("Aureo_data2.csv"),  file("data/deseq_fonc3.R") )
}
