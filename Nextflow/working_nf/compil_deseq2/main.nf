nextflow.enable.dsl=2

process DOWNLOAD_SRA {

    tag "$srr_id"
    container 'pierrejeangouze/sra_toolkit'
    errorStrategy 'retry'
    maxRetries 2

    input:
    val srr_id

    output:
    path "*.fastq.gz"

    script:
    """
    # Telechargement
    fasterq-dump --threads 12 --split-files --progress ${srr_id}
    # Compression
    gzip *.fastq
    """
}


process TRIM_ALL_READS {
    tag "${fastq_file.simpleName}"
    container "bagranier/trim_galore_cutadapt"

    input:
    path fastq_file

    output:
    path "*.fq.gz"

    script: 
    """
    trim_galore -q 20 --phred33 --length 25 ${fastq_file}
    """
}

//--------------------------------------------------------------
// 1. DOWNLOAD REFERENCE
//--------------------------------------------------------------
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



//--------------------------------------------------------------
// 2. BOWTIE INDEX CREATION (all index files in a folder)
//--------------------------------------------------------------
process INDEX_CREATION {
    tag "index"
    container 'vmichelet/bowtie:0.12.7-samtools'
    containerOptions '--entrypoint ""'

    input:
    path fasta_file

    output:
    path "ref_index", emit: index_dir

    script:
    """
    mkdir ref_index
    bowtie-build ${fasta_file} ref_index/ref_index
    """
}



//--------------------------------------------------------------
// 3. ALIGN READS USING PROCESS SUBSTITUTION + PIPED SORT
//--------------------------------------------------------------
process ALIGN {
    tag "align"
    container 'vmichelet/bowtie:0.12.7-samtools'
    containerOptions '--entrypoint ""'
    cpus 4

    input:
    path fastq_file
    path index_dir

    output:
    path "*.bam"

    script:
    """
    bowtie -p ${task.cpus} -S ${index_dir}/ref_index <(gunzip -c ${fastq_file}) \
        | samtools sort -@ ${task.cpus} -o ${fastq_file.simpleName}.bam

    samtools index ${fastq_file.simpleName}.bam
    """
}

process FEATURECOUNTS {

    container 'pierrejeangouze/featurecounts:1.4.6-p3'

    input:
        path bam_list
        path gff

    output:
        path "counts.txt"

    script:
    """
    featureCounts \
        -t gene \
        -g ID \
        -F GFF \
        -T 4 \
        -a $gff \
        -o counts.txt ${bam_list.join(' ')}
    """
}

process deseq2 {
    tag "Running"
    container 'vmichelet/r341_desq2'

    input:
    path counts
    path reference
    path tsv

    output:
    path "*_simple.pdf"
    path "*_annotated.pdf"
    path "DESeq2_results_clean.csv"

    script:
    """
    Rscript - << 'EOF'
    ############################################################
    ## BLOC 1 — Import counts.txt + préparation DESeq2
    ############################################################

    library(DESeq2)
    library(ggplot2)

    ## 1. Lecture du fichier featureCounts
    raw_counts <- read.table(
      "${counts}",
      header = TRUE,
      row.names = 1,
      sep = "\t",
      check.names = FALSE
    )

    ## 2. Sélection des colonnes d’échantillons
    ## IMPORTANT : chez toi, elles sont bien colonnes 6 à 11
    counts <- raw_counts[, 6:11]

    ## 3. Vérification
    print(colnames(counts))

    ## 4. Conversion en matrice entière
    counts <- as.matrix(counts)
    storage.mode(counts) <- "integer"

    ## 5. Récupération des noms d’échantillons
    samples <- colnames(counts)

    ## 6. Définition des conditions (ordre = celui des colonnes)
    condition <- c(
      "control","control","control",
      "persister","persister","persister"
    )

    ## 7. Construction de colData
    colData <- data.frame(
      row.names = samples,
      condition = factor(condition)
    )

    ## 8. Création de l’objet DESeq2
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = colData,
      design = ~ condition
    )

    ## 9. Filtrage minimal (gènes exprimés)
    dds <- dds[rowSums(counts(dds)) > 1, ]

    ## 10. Lancement de DESeq2
    dds <- DESeq(dds)

    ## 11. Extraction des résultats
    res <- results(dds)
    df <- as.data.frame(res)

    ############################################################
    ## BLOC 1 TERMINÉ — df contient vos résultats DESeq2
    ############################################################

    ############################################################
    ## BLOC 2 — Import du GFF + extraction des locus tags
    ############################################################

    ## 1. Lire le fichier GFF
    gff <- read.table(
      "${reference}",
      sep = "\t",
      header = FALSE,
      comment.char = "#",
      stringsAsFactors = FALSE
    )

    ## 2. Extraire la colonne des attributs (colonne 9)
    attributes <- gff\$V9

    ## 3. Fonction robuste pour extraire ID=SAUSA300_xxxx
    extract_locus <- function(x) {
      sub(".*ID=([^;]+).*", "\\1", x)
    }

    ## 4. Extraction effective
    gff_ids <- extract_locus(attributes)

    ## 5. Nettoyage : si l’extraction ratée → NA
    gff_ids[gff_ids == attributes] <- NA

    ## 6. Facultatif : data.frame structuré pour le génome
    gff_df <- data.frame(
      locus_tag = gff_ids,
      chr = gff\$V1,
      start = gff\$V4,
      end = gff\$V5,
      strand = gff\$V7,
      stringsAsFactors = FALSE
    )

    ## 7. Vérification
    head(gff_df)

    ############################################################
    ## BLOC 2 TERMINÉ — locus_tag = SAUSA300_xxxx récupérés
    ############################################################


    ############################################################
    ## BLOC 3 — MAPPING AUREOWIKI (final et corrigé)
    ############################################################

    ## 1. Lecture du fichier AureoWiki COL
    mapping <- read.delim(
      "${tsv}",
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    )

    ## 2. Vérification des colonnes
    colnames(mapping)

    ## 3. Sélection des colonnes utiles :
    ##   - locus.tag  : identifiant de gène (équivalent SAOUHSC)
    ##   - symbol     : nom de gène
    mapping_clean <- mapping[, c("locus.tag", "symbol")]

    ## 4. Renommage propre pour la suite du workflow
    colnames(mapping_clean) <- c("Locus_Tag", "Gene_Name")

    ## 5. Vérification finale
    head(mapping_clean)

    ############################################################
    ## → mapping_clean\$Locus_Tag : SACOLxxxx
    ## → mapping_clean\$Gene_Name : dnaA, dnaN, recF, …
    ############################################################

    ############################################################
    ## BLOC 4 — Liste locale des gènes de la traduction
    ## (compatible R 3.4.1, aucune dépendance KEGG)
    ############################################################

    # Gènes codant les protéines du ribosome 30S
    ribosomal_30S <- c(
      "rpsA","rpsB","rpsC","rpsD","rpsE","rpsF","rpsG",
      "rpsH","rpsI","rpsJ","rpsK","rpsL","rpsM","rpsN",
      "rpsO","rpsP","rpsQ","rpsR","rpsS"
    )

    # Gènes codant les protéines du ribosome 50S
    ribosomal_50S <- c(
      "rplA","rplB","rplC","rplD","rplE","rplF","rplG",
      "rplH","rplI","rplJ","rplK","rplL","rplM","rplN",
      "rplO","rplP","rplQ","rplR","rplS","rplT","rplU",
      "rplV","rplW","rplX","rplY"
    )

    # Facteurs d'initiation / élongation de la traduction
    translation_factors <- c("tsf","tuf","fusA","infA","infB","infC")

    # Liste finale
    ribosomal_gene_names <- c(
      ribosomal_30S,
      ribosomal_50S,
      translation_factors
    )

    # Vérification
    ribosomal_gene_names
    length(ribosomal_gene_names)

    ############################################################
    ## → Environ 55–60 gènes
    ## → Utilisables directement pour les MA-plots
    ############################################################

    ############################################################
    ## BLOC 5 — Construction du tableau final avec DESeq2 + ribosomal
    ############################################################

    ## 1. Convertir les résultats DESeq2 en data.frame
    res_df <- as.data.frame(res)
    res_df\$geneID <- rownames(res_df)

    ## 2. Nettoyer les IDs (supprimer le préfixe "gene-")
    res_df\$Locus_Tag <- sub("^gene-", "", res_df\$geneID)

    ## 3. Liste des gènes ribosomiques SAOUHSC (USA300/COL cluster)
    ribosomal_saouhsc <- c(
      "SAOUHSC_02120","SAOUHSC_02121","SAOUHSC_02122","SAOUHSC_02123",
      "SAOUHSC_02124","SAOUHSC_02125","SAOUHSC_02126","SAOUHSC_02127",
      "SAOUHSC_02128","SAOUHSC_02129","SAOUHSC_02130","SAOUHSC_02131",
      "SAOUHSC_02132","SAOUHSC_02133","SAOUHSC_02134","SAOUHSC_02135",
      "SAOUHSC_02136","SAOUHSC_02137","SAOUHSC_02138","SAOUHSC_02139",
      "SAOUHSC_02140","SAOUHSC_02141","SAOUHSC_02142","SAOUHSC_02143",
      "SAOUHSC_02144","SAOUHSC_02145","SAOUHSC_02146","SAOUHSC_02147",
      "SAOUHSC_02148","SAOUHSC_02149","SAOUHSC_02150"
    )

    ## 4. Annotation ribosomique
    res_df\$is_ribosomal <- res_df\$Locus_Tag %in% ribosomal_saouhsc

    ## 5. Vérification
    table(res_df\$is_ribosomal)

    ## 6. Export fichier annoté
    write.csv(
      res_df,
      file = "DESeq2_results_clean.csv",
      row.names = FALSE
    )

    ############################################################
    ## BLOC 5 TERMINÉ — res_df est prêt pour les figures du TP
    ############################################################

    ############################################################
    ## BLOC 6 — Génération des figures du TP (simple + annotée)
    ############################################################

    library(ggplot2)

    ############################################################
    ## 6.1 — MA plot simple (comme dans DESeq2)
    ############################################################

    pdf("MA_plot_simple.pdf", width = 7, height = 6)

    plotMA(
      res,
      ylim = c(-4,4),
      alpha = 0.05,
      main = "MA plot - Persister vs Control"
    )

    dev.off()


    ############################################################
    ## 6.2 — MA plot annoté (gènes ribosomiques en rouge)
    ############################################################

    pdf("MA_plot_annotated.pdf", width = 7, height = 6)

    ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
      geom_point(
        aes(color = is_ribosomal),
        alpha = 0.6,
        size = 1.5
      ) +
      scale_x_log10() +
      scale_color_manual(
        values = c("FALSE" = "grey60", "TRUE" = "red"),
        labels = c("Other genes", "Ribosomal genes")
      ) +
      theme_minimal(base_size = 12) +
      labs(
        title = "MA plot with ribosomal gene annotation",
        x = "Mean expression (log10)",
        y = "log2 Fold Change"
      )

    dev.off()

    ############################################################
    ## BLOC 6 TERMINÉ — Les PDF sont générés dans ton répertoire
    ############################################################

    EOF
    """
}


workflow {
    // Liste definie directement ici
    srr_list = [
        'SRR10379721',
        'SRR10379722',
        'SRR10379723',
        'SRR10379724',
        'SRR10379725',
        'SRR10379726'
    ]

    // Afficher la liste utilisee
    println "Liste des echantillons : ${srr_list.join(', ')}"

    // Creer un canal qui emet chaque SRR separement
    srr_ch = Channel.from(srr_list)

    // Lancer une tache pour chaque SRR
    fastq_files = DOWNLOAD_SRA(srr_ch)
    (fasta, gff) = DOWNLOAD_REFERENCE()
    index = INDEX_CREATION(fasta)
    fastq_ch = TRIM_ALL_READS(fastq_files)
    bam = ALIGN(fastq_ch, index)
    counts = FEATURECOUNTS(bam.collect(), gff)
    deseq2(counts,gff, file("GeneSpecificInformation_COL.tsv") )
}
