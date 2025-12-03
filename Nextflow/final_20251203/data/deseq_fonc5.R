
args <- commandArgs(trailingOnly = TRUE)
counts  <- args[1]   # counts.txt
aureo  <- args[2]   # liste gèenes
res_exp <- args[3] # résultats des chercheurs

#############################################
### 
#############################################

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
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"


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

# -- 1. Charger Aureo_data.csv (séparateur ;) --
raw <- read.csv2(aureo, header = TRUE, stringsAsFactors = FALSE)

# Le fichier contient déjà les colonnes :
#   GeneID      -> ex : "SAOUHSC_00017"
#   GeneName    -> ex : "rplI"

translation_genes <- raw

# -- 2. Normaliser les GeneID pour correspondre à DESeq2 --
translation_genes$GeneID <- toupper(trimws(translation_genes$GeneID))

# -- 3. Récupérer les résultats DESeq2 --
res_df <- as.data.frame(res)

# Ajouter colonne GeneID nettoyée
res_df$GeneID <- toupper(gsub("^gene-", "", rownames(res_df)))

# -- 4. Filtrer pour ne garder que les gènes présents dans Aureo_data --
sub <- merge(res_df, translation_genes, by = "GeneID")


# -----------------------------
# 1. Liste des gènes à mettre en avant
# -----------------------------
highlight_genes <- c("frr","infA","tsf","infC","infB","pth")


# Merge final propre
sub <- merge(res_df, translation_genes, by = "GeneID")

# -----------------------------
# 3. Filtrage gènes de traduction
# -----------------------------
sub <- res_df[res_df$GeneID %in% translation_genes$GeneID, ]
sub <- merge(sub, translation_genes, by = "GeneID")

# -----------------------------
# 4. Définir AA-tRNA synthetases
# -----------------------------
aaRS_list <- c("ileS","leuS","metG","pheS","pheT","proS",
               "serS","thrS","trpS","tyrS","valS","alaS",
               "argS","asnS","aspS","cysS","glnS","glyS","hisS","lysS")
sub$aaRS <- sub$GeneName %in% aaRS_list

# -----------------------------
# 5. Définir couleur selon significativité (padj < 0.05)
# -----------------------------
sub$color <- ifelse(sub$padj < 0.05, "red", "grey")

# -----------------------------
# 6. Ligne basale
# -----------------------------
basal_line <- 0

# # -----------------------------
# # 7. MA-plot complet avec traits et labels
# # -----------------------------
plot_MA_translation <- function(sub, fig_title = "MA-plot") {

#png("MA_translation_genes_final.png", width = 1400, height = 1000)
png(paste0(gsub(" ", "_", fig_title), ".png"), width = 700, height = 700)  # carré
#pdf("MA_plot_translationgenes_FINAL.pdf", width = 10, height = 10)
  # Dataframe
  plot_df <- data.frame(
    log2BaseMean     = log2(sub$baseMean),
    log2FoldChange   = sub$log2FoldChange,
    signif           = ifelse(sub$color == "red", "Significant", "Non-Significant"),
    is_aaRS          = sub$aaRS,
    GeneName         = sub$GeneName
  )

  highlight_df <- subset(plot_df, GeneName %in% highlight_genes)

  # ----- BASE GGplot -----
p <- ggplot(plot_df, aes(log2BaseMean, log2FoldChange)) +

  # -------------------------------------------------------------------------
  # Points gris et rouges
  geom_point(aes(color = signif), size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c("Non-Significant" = "grey70",
               "Significant"     = "red"),
    name = "Genes"
  ) +

  # AA-tRNA synthétases
  geom_point(
    data = subset(plot_df, is_aaRS),
    color = "red", size = 3.5
  ) +
  geom_point(
    data = subset(plot_df, is_aaRS),
    aes(shape = "AA_tRNA_synthetases"),
    size = 4, stroke = 2.2, color = "black"
  ) +
  scale_shape_manual(values = c("AA_tRNA_synthetases" = 1), name = "") +

  # Segments et labels
  geom_segment(
    data = highlight_df,
    aes(xend = log2BaseMean + 0.8, yend = log2FoldChange + 0.8),
    size = 0.9
  ) +
  geom_text(
    data = highlight_df,
    aes(x = log2BaseMean + 0.9, y = log2FoldChange + 0.9, label = GeneName),
    fontface = "italic", size = 5
  ) +

  geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) +

  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-6, 5),
    breaks = seq(-6, 5, by = 1),
    expand = c(0, 0)
  ) +

  # -------------------------------------------------------------------------
  # carré parfait + légende dans le cadre
  coord_fixed(ratio = 2) +
  theme_bw() +
  ggtitle(fig_title)+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),

    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),

    legend.position = c(0.05, 0.1),      # <<< position interne (en haut à gauche)
    legend.justification = c("left"),
    legend.box = "horizontal",
    legend.background = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),

    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_text(
        size = 22,       # << taille
        #face = "bold",   # optionnel, pour un titre plus propre
        hjust = 0.5      # << CENTRÉ
    ),
  )
  print(p)

dev.off()
}



plot_all_genes <- function(res, fig_title = "MA-plot_all") {

png(paste0(gsub(" ", "_", fig_title), ".png"), width=1000, height=800)
#png("okok.png", width=1000, height=800)

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
print(
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
    ggtitle(paste0("MA-plot with triangles for out-of-range values _", fig_title)) +
    xlab("baseMean (log10 scale)") +
    ylab("log2 Fold Change")+
    theme(
        plot.title = element_text(
            size = 20,      # ← augmente la taille
            hjust = 0.5     # ← CENTRE le titre
        )
    )
)
dev.off()
}

#############################################
###       VOLCANO PLOT                    ###
#############################################

plot_volcano <- function(res, fig_title = "Volcano_plot") {

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

png(paste0(gsub(" ", "_", fig_title), ".png"), width = 1200, height = 1000)
print(
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
        plot.title = element_text(
            size = 20,      # ← augmente la taille
            hjust = 0.5     # ← CENTRE le titre
        )
    ) +
    
    labs(
        title = "Volcano Plot (DESeq2)",
        x = "log2 Fold Change (-5,5)",
        y = "-log10(p-adj) (-1,75)",
        color = "Gene Group"
    )
)
dev.off()
}

######

plot_MA_translation(sub, "MA-plot_repro") 

plot_all_genes(res, "MA_plot_all_genes_repro")

plot_volcano(res, "Volcano_repro")

#### Figures de l'expérience

# plot_MA_translation(res_exp, "MA-plot_repro") 

# plot_all_genes(res_exp, "MA_plot_all_genes_repro")

# plot_volcano(res_exp, "Volcano_repro")

write.csv2(sub, "Resultats_deseq2.csv", row.names = FALSE)
