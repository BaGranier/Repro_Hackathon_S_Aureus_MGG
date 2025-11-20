nextflow.enable.dsl=2

process DESeq2 {

    tag "Running-DESeq2"
    container 'vmichelet/r341_desq2'

    input:
    path counts_file
    path gff_file
    path aureo_file

    output:
    path "MA_plot_simple.pdf"
    path "MA_plot_annotated.pdf"
    path "MA_plot_ggplot2.pdf"
    path "DESeq2_results_clean.csv"

    script:
    """
    mkdir -p data
    cp ${counts_file} data/counts.txt
    cp ${gff_file} data/reference.gff
    cp ${aureo_file} data/GeneSpecificInformation_COL.tsv

    Rscript - << 'EOF'

    #############################
    ## BLOC 1 : Import DESeq2
    #############################

    library(DESeq2)
    library(ggplot2)
    library(dplyr)

    raw_counts <- read.table(
      "data/counts.txt",
      header = TRUE,
      row.names = 1,
      sep = "\\t",
      check.names = FALSE
    )

    counts <- raw_counts[, 6:11]
    counts <- as.matrix(counts)
    storage.mode(counts) <- "integer"

    samples <- colnames(counts)
    condition <- c("control","control","control",
                   "persister","persister","persister")

    colData <- data.frame(
      row.names = samples,
      condition = factor(condition)
    )

    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = colData,
      design = ~ condition
    )

    dds <- dds[rowSums(counts(dds)) > 1, ]
    dds <- DESeq(dds)
    res <- results(dds)
    df <- as.data.frame(res)


    #############################
    ## BLOC 2 : Import GFF
    #############################

    gff <- read.table(
      "data/reference.gff",
      sep = "\\t",
      header = FALSE,
      comment.char = "#",
      stringsAsFactors = FALSE
    )
    attributes <- gff$V9
    extract_locus <- function(x) sub(".*ID=([^;]+).*", "\\\\1", x)
    gff_ids <- extract_locus(attributes)
    gff_ids[gff_ids == attributes] <- NA


    #############################
    ## BLOC 3 : AureoWiki mapping
    #############################

    mapping <- read.delim(
      "data/GeneSpecificInformation_COL.tsv",
      header = TRUE,
      sep = "\\t",
      fill = TRUE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    )

    mapping_clean <- mapping[, c("locus.tag", "symbol")]
    colnames(mapping_clean) <- c("Locus_Tag", "Gene_Name")


    #############################
    ## BLOC 4 : Ribosomal genes
    #############################

    ribosomal_30S <- c("rpsA","rpsB","rpsC","rpsD","rpsE","rpsF","rpsG",
                       "rpsH","rpsI","rpsJ","rpsK","rpsL","rpsM","rpsN",
                       "rpsO","rpsP","rpsQ","rpsR","rpsS")

    ribosomal_50S <- c("rplA","rplB","rplC","rplD","rplE","rplF","rplG",
                       "rplH","rplI","rplJ","rplK","rplL","rplM","rplN",
                       "rplO","rplP","rplQ","rplR","rplS","rplT","rplU",
                       "rplV","rplW","rplX","rplY")

    translation_factors <- c("tsf","tuf","fusA","infA","infB","infC")

    ribosomal_gene_names <- c(ribosomal_30S, ribosomal_50S, translation_factors)


    #############################
    ## BLOC 5 : Annot final
    #############################

    df$geneID <- rownames(df)
    df$Locus_Tag <- sub("^gene-", "", df$geneID)

    df <- df %>%
        left_join(mapping_clean, by = "Locus_Tag")

    df$is_ribo <- df$Gene_Name %in% ribosomal_gene_names

    write.csv(df, "DESeq2_results_clean.csv", row.names = FALSE)


    #############################
    ## BLOC 6 : Plots PDF
    #############################

    ## MA plot simple
    pdf("MA_plot_simple.pdf", width = 7, height = 6)
    plotMA(res, ylim = c(-4,4), alpha = 0.05)
    dev.off()

    ## MA plot annoté
    pdf("MA_plot_annotated.pdf", width = 7, height = 6)
    ggplot(df, aes(x = baseMean, y = log2FoldChange)) +
      geom_point(aes(color = is_ribo), alpha = 0.6, size = 1.5) +
      scale_x_log10() +
      scale_color_manual(values = c("grey60","red")) +
      theme_minimal()
    dev.off()

    ## MA plot ggplot2 avec 10^n
    df$baseMean_log10 <- log10(df$baseMean + 1)
    pdf("MA_plot_ggplot2.pdf", width = 7, height = 6)
    ggplot(df, aes(x = baseMean_log10, y = log2FoldChange)) +
      geom_point(aes(color = is_ribo), alpha = 0.6) +
      scale_x_continuous(
        breaks = 0:6,
        labels = scales::math_format(10^.x)
      ) +
      theme_bw()
    dev.off()

    EOF
    """
}

workflow {
    DESeq2(
        file("counts.txt"),
        file("reference.gff"),
        file("GeneSpecificInformation_COL.tsv")
    )
}