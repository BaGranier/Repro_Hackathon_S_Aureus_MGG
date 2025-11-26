nextflow.enable.dsl=2

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
    deseq2(file("counts.txt"),file("reference.gff"), file("data/GeneSpecificInformation_COL.tsv"), file("data/deseq.R") )
}
