#!/usr/bin/env nextflow

process FEATURECOUNTS {

    container 'pierrejeangouze/featurecounts:1.4.6-p3'

    input:
    path bam
    path gff

    output:
    path "counts.txt"

    """
    featureCounts \
        -t gene \
        -g ID \
        -F GFF \
        -T 4 \
        -a $gff \
        -o counts.txt $bam
    """
}

workflow {

    // Channel avec ton(s) fichier(s) BAM
    bam_ch = Channel.fromPath("SRR*.bam")

    // Channel avec le fichier GFF
    gff_ch = Channel.fromPath("reference.gff")

    FEATURECOUNTS(bam_ch, gff_ch)
}

