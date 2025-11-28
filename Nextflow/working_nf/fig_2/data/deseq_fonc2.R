

args <- commandArgs(trailingOnly = TRUE)
counts  <- args[1]   # counts.txt
aureo  <- args[2]   # counts.txt

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


#############################################
### VOLCANO PLOT — compatible R 3.4.1     ###
### ggplot2 2.2.1 — SANS dplyr, SANS ggrepel
#############################################
library(ggplot2)

df <- as.data.frame(res)

# Retirer NA
df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), ]

# Seuils
padj_threshold  <- 0.05
logfc_threshold <- 1

# Catégorie Up / Down / NS
df$signif <- ifelse(df$padj < padj_threshold & df$log2FoldChange >  logfc_threshold, "Up",
             ifelse(df$padj < padj_threshold & df$log2FoldChange < -logfc_threshold, "Down",
             "NS"))

# Sélection des gènes à annoter (top 10)
df_sorted <- df[order(df$padj), ]
df_labels <- df_sorted[df_sorted$signif != "NS", ]
df_labels <- df_labels[1:min(10, nrow(df_labels)), ]


png("Volcano_ggplot_limited_R341.png", width = 1200, height = 1000)

ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    
    # Points
    geom_point(aes(color = signif), alpha = 0.7, size = 1.8) +
    
    # Couleurs manuelles compatibles ggplot2 2.2.1
    scale_color_manual(values = c(
        "NS"   = "grey60",
        "Up"   = "firebrick2",
        "Down" = "dodgerblue2"
    )) +
    
    # Lignes de seuil
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
               color = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_threshold),
               color = "black", linetype = "dashed") +
    
    # Labels simples (sans ggrepel)
    geom_text(data = df_labels,
              aes(label = rownames(df_labels)),
              size = 3,
              vjust = -0.5) +

    # 🟢 LIMITES D’AXES (NEW)
    #coord_cartesian(xlim = c(-5, 5), ylim = c(-1, 75)) +
    
    theme_bw(base_size = 16) +
    theme(
        plot.title = element_text(face = "bold", size = 20),
        legend.position = "right"
    ) +
    
    labs(
        title = "Volcano Plot (DESeq2)",
        x = "log2 Fold Change (-5,5)",
        y = "-log10(p-adj) (-1,75)",
        color = "Gene Group"
    )

dev.off()


# # Charger le fichier Aureo_data.csv contenant tous les gènes de traduction
# raw <- read.csv2(aureo, header = TRUE, stringsAsFactors = FALSE)

# # La colonne s'appelle GeneID.GeneName, on la sépare en GeneID et GeneName
# translation_genes <- raw[, c("GeneID", "GeneName")]


# # Normaliser les GeneID pour correspondre à res$GeneID
# translation_genes$GeneID <- toupper(trimws(translation_genes$GeneID))

# # Liste des gènes à mettre en avant
# highlight_genes <- c("frr","infA","tsf","infC","infB","pth")

# # Nettoyage GeneID
# res$GeneID <- gsub("^gene-", "", rownames(res))
# res_df <- as.data.frame(res)
# translation_genes$GeneID <- toupper(trimws(translation_genes$GeneID))


# # 1. Charger Aureo
# raw <- read.csv2(aureo, header = TRUE, stringsAsFactors = FALSE)
# translation_genes <- raw[, c("GeneID", "GeneName")]

# # 2. Rendre les GeneID comparables
# translation_genes$GeneID <- toupper(trimws(translation_genes$GeneID))

# # 3. Convertir les résultats DESeq2 en data.frame
# res_df <- as.data.frame(res)

# # 4. Sauvegarder les ID bruts
# res_df$OriginalID <- rownames(res_df)

# # -----------------------------
# # 1. Liste des gènes à mettre en avant
# # -----------------------------
# highlight_genes <- c("frr","infA","tsf","infC","infB","pth")

# # -----------------------------
# # 2. Nettoyage GeneID
# # -----------------------------
# res$GeneID <- gsub("^gene-", "", rownames(res))
# res_df <- as.data.frame(res)
# translation_genes$GeneID <- toupper(trimws(translation_genes$GeneID))

# # -----------------------------
# # 3. Filtrage gènes de traduction
# # -----------------------------
# sub <- res_df[res_df$GeneID %in% translation_genes$GeneID, ]
# sub <- merge(sub, translation_genes, by = "GeneID")

# # -----------------------------
# # 4. Définir AA-tRNA synthetases
# # -----------------------------
# aaRS_list <- c("ileS","leuS","metG","pheS","pheT","proS",
#                "serS","thrS","trpS","tyrS","valS","alaS",
#                "argS","asnS","aspS","cysS","glnS","glyS","hisS","lysS")
# sub$aaRS <- sub$GeneName %in% aaRS_list

# # -----------------------------
# # 5. Définir couleur selon significativité (padj < 0.05)
# # -----------------------------
# sub$color <- ifelse(sub$padj < 0.05, "red", "grey")

# # -----------------------------
# # 6. Ligne basale
# # -----------------------------
# basal_line <- 0

# # -----------------------------
# # 7. MA-plot complet avec traits et labels
# # -----------------------------
# png("MA_translation_genes_final.png", width = 1400, height = 1000)
# plot(
#   log2(sub$baseMean),
#   sub$log2FoldChange,
#   pch = 16,
#   cex = 2,               # taille des points
#   col = sub$color,
#   xlab = "log2(Base Mean)",
#   ylab = "log2 Fold Change",
#   main = "MA-plot – Translation Genes",
#   ylim = c(-6, 5),
#   xlim = c(0, 20),
#   xaxt = "n",
#   yaxt = "n"
# )

# # Ligne basale
# abline(h = basal_line, lty = 2, col = "black")

# # Axes personnalisés
# axis(side = 1, at = seq(0, 20, by = 2))
# axis(side = 2, at = seq(-6, 5, by = 1))

# # Ajouter contours noirs pour AA-tRNA synthetases
# points(
#   log2(sub$baseMean[sub$aaRS]),
#   sub$log2FoldChange[sub$aaRS],
#   pch = 1,   # cercle vide
#   cex = 2.5,
#   lwd = 2,
#   col = "black"
# )

# # Ajouter traits verticaux pour les gènes à mettre en avant
# highlight_points <- sub[sub$GeneName %in% highlight_genes, ]
# for(i in 1:nrow(highlight_points)){
#   segments(
#     x0 = log2(highlight_points$baseMean[i]),
#     y0 = highlight_points$log2FoldChange[i],
#     x1 = log2(highlight_points$baseMean[i]),
#     y1 = highlight_points$log2FoldChange[i] + 1,  # longueur du trait
#     col = "red",
#     lwd = 1.5,
#     lty = 2
#   )
# }

# # Ajouter labels au bout des traits
# text(
#   log2(highlight_points$baseMean + 1),
#   highlight_points$log2FoldChange + 1.1,
#   labels = highlight_points$GeneName,
#   pos = 3,
#   cex = 1,
#   col = "red"
# )

# # Légende
# legend(
#   "topright",
#   legend = c("Significant", "Non-significant", "aaRS", "Highlighted genes"),
#   col = c("red", "grey", "black", "red"),
#   pch = c(16,16,1,NA),
#   lty = c(NA,NA,NA,2),
#   pt.cex = c(2,2,2.5,NA),
#   pt.lwd = c(1,1,2,1)
# )

# dev.off()