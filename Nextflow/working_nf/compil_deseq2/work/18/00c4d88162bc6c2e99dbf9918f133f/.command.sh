#!/bin/bash -ue
featureCounts         -t gene         -g ID         -F GFF         -T 4         -a reference.gff         -o counts.txt SRR10379726_1_trimmed.bam SRR10379723_1_trimmed.bam SRR10379722_1_trimmed.bam SRR10379725_1_trimmed.bam SRR10379724_1_trimmed.bam SRR10379721_1_trimmed.bam
