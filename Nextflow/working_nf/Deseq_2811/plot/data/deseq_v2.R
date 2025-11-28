

args <- commandArgs(trailingOnly = TRUE)
counts  <- args[1]   # counts.txt

#1 
library(DESeq2)
library(ggplot2)
# Lire le counts.txt
counts_raw <- read.table(counts, header=TRUE, sep="\t",
                         comment.char="#", stringsAsFactors=FALSE, check.names=FALSE)

desired_order <- c(
  "SRR10379721_1_trimmed.bam",  # IP1
  "SRR10379722_1_trimmed.bam",  # IP2
  "SRR10379723_1_trimmed.bam",  # IP3
  "SRR10379724_1_trimmed.bam",  # Ctrl1
  "SRR10379725_1_trimmed.bam",  # Ctrl2
  "SRR10379726_1_trimmed.bam"   # Ctrl3
)


# Colonnes des counts
start_col <- 7
count_matrix <- counts_raw[, start_col:ncol(counts_raw)]
count_matrix <- count_matrix[, desired_order]

# Gene IDs comme noms de lignes
rownames(count_matrix) <- counts_raw$Geneid

# Convertir en matrice entière
count_matrix <- as.matrix(apply(count_matrix, 2, as.integer))

# Créer coldata
sample_names <- colnames(count_matrix)
coldata <- data.frame(
  row.names = sample_names,
  condition = factor(c("persister","persister","persister",
                       "control","control","control"))
)


# Créer l'objet DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ condition)

# Filtrer les gènes avec peu de counts (au total < 10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


#3
# Lancer l'analyse DESeq2
dds <- DESeq(dds)

# Obtenir les résultats (padj < 0.05)
res <- results(dds, alpha=0.05)

# Trier les gènes par p-value ajustée
resOrdered <- res[order(res$padj),]

#4
# MA-plot de tous les gènes
# png("MAplot_all_genes.png", width=800, height=600)
# plotMA(res, ylim=c(-4,4), main="MA-plot (all genes)")
# plot.window(xlim = c(0, 1e6), ylim = c(-4, 4),  log = "x")
# dev.off()



png("MAplot_custom.png", width=1000, height=800)

df <- as.data.frame(res)
df <- df[df$baseMean <= 1e6, ]

# ---- Config ----
ylim_min <- -4
ylim_max <- 4

# couleur DESeq2-like
df$col <- ifelse(df$padj < 0.05, "red", "black")

# ---- Séparer les 3 catégories ----
df_inside  <- df[df$log2FoldChange <= ylim_max & df$log2FoldChange >= ylim_min, ]
df_top     <- df[df$log2FoldChange >  ylim_max, ]
df_bottom  <- df[df$log2FoldChange <  ylim_min, ]

# ---- Plot ----
ggplot() +
    # Points "normaux"
    geom_point(data = df_inside,
               aes(x = baseMean, y = log2FoldChange, color = col),
               alpha = 0.6, size = 1) +

    # Triangles du haut
    geom_point(data = df_top,
               aes(x = baseMean, y = ylim_max, color = col),
               shape = 24, size = 2, fill = "white") +     # triangle haut

    # Triangles du bas
    geom_point(data = df_bottom,
               aes(x = baseMean, y = ylim_min, color = col),
               shape = 25, size = 2, fill = "white") +     # triangle bas

    scale_color_identity() +
    scale_x_log10(limits = c(1, 1e6)) +
    ylim(ylim_min, ylim_max) +
    ggtitle("MA-plot with triangles for out-of-range values") +
    xlab("baseMean (log10 scale)") +
    ylab("log2 Fold Change")

dev.off()

############################################################
## 6.2 — MA plot annoté (gènes ribosomiques en rouge)
############################################################

if (!requireNamespace("EnrichmentBrowser", quietly = TRUE)) {
    BiocManager::install("EnrichmentBrowser")
}
library(EnrichmentBrowser)
# Obtenir la liste des gènes ribosomiques pour l'espèce donnée
library(EnrichmentBrowser)

list_kegg <- list("sao03010","sao00970")
all_gene <-  list("ID" = list(), "NAME" = list())

for (KEgg in list_kegg) {
  gene_keg <- KEGGREST::keggGet(KEgg)[[1]]
  print(KEgg)
  
  for (i in seq(1,length((gene_keg$GENE)), 2)) {
    
    all_gene$ID = append(all_gene$ID, gene_keg$GENE[i])
    all_gene$NAME = append(all_gene$NAME, gene_keg$GENE[i+1])
    
  }
  
}

# Missing gene

all_gene$ID = append(all_gene$ID, "SAOUHSC_01203")
all_gene$NAME = append(all_gene$NAME, "rnc; ribonuclease III [KO:K03685] [EC:3.1.26.3]")

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