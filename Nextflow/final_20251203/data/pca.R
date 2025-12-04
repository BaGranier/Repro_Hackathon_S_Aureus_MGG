
args <- commandArgs(trailingOnly = TRUE)
counts  <- args[1]   # counts.txt


library(DESeq2)
library(ggplot2)
library(ggrepel)

# -----------------------------
# 1. Importer counts.txt
# -----------------------------
counts <- read.table(counts, header = TRUE, row.names = 1, sep = "\t")

desired_order <- c(
  "SRR10379721_1_trimmed.bam",  # IP1
  "SRR10379722_1_trimmed.bam",  # IP2
  "SRR10379723_1_trimmed.bam",  # IP3
  "SRR10379724_1_trimmed.bam",  # Ctrl1
  "SRR10379725_1_trimmed.bam",  # Ctrl2
  "SRR10379726_1_trimmed.bam"   # Ctrl3
)

# Sélection des colonnes d'échantillons
sample_cols <- grep("bam$", colnames(counts))
count_matrix <- counts[, sample_cols]
count_matrix <- count_matrix[, desired_order]

# -----------------------------
# 2. Créer le metadata
# -----------------------------
metadata <- data.frame(
  Sample = colnames(count_matrix),
  Condition = c("persister","persister","persister","control","control","control")
)
rownames(metadata) <- metadata$Sample
metadata$Condition <- factor(metadata$Condition, levels = c("control","persister"))

# -----------------------------
# 3. Création DESeq2 + rlog
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ Condition)
dds <- DESeq(dds)
rld <- rlog(dds)

# -----------------------------
# 4. Extraction des données PCA
# -----------------------------
pca_data <- plotPCA(rld, intgroup = "Condition", returnData = TRUE)
percent_variance <- round(100 * attr(pca_data, "percentVar"))

# Couper les noms avant le "."
pca_data$name <- sapply(pca_data$name, function(x) strsplit(x, "\\.")[[1]][1])

# -----------------------------
# 5. Création du plot ggplot2 avec ggrepel
# -----------------------------
png("pca_repro.png", width = 700, height = 700)  # carré
ACP <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = name)) +
        geom_point(size = 4, alpha = 0.8) +
        geom_text_repel(size = 3, max.overlaps = Inf) +
        labs(title = "Principal Component Analysis of samples (PCA)",
             x = paste0("PC1: ", percent_variance[1], "% variance"),
             y = paste0("PC2: ", percent_variance[2], "% variance"),
             color = "Condition") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom",
              plot.title = element_text(size = 20, hjust = 0.5),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16))

print(ACP)
dev.off()
# -----------------------------
# 6. Sauvegarde du plot
# -----------------------------
# ggsave("ACP_samples.png", plot = ACP, width = 10, height = 8, dpi = 300)

# cat("Plot ACP généré : ACP_samples.png\n")