nextflow.enable.dsl=2

process deseq2 {
    tag "Running"
    container 'vmichelet/r341_desq2'

    input:
    path counts
    path aureo
    path res_exp
    file deseq_script

    output:
    path "*.png"

    script:
    """
    Rscript ${deseq_script} ${counts} ${aureo}
    """
}

workflow {
    deseq2(file("counts.txt"), file("data/Aureo_data2.csv"), file("data/GSE139659_IPvsctrl.complete.csv"),  file("data/deseq_fonc4.R") )
}
