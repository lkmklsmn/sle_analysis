# Load R libs ####
library(GEOquery)
library(ggplot2)
library(limma)


# Download data from GEO ####
chiche <- getGEO("GSE49454", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
banchereau <- getGEO("GSE65391", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
#hong <- getGEO("GSE108497", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
# Hong doesnt download correctly
petri_1 <- getGEO("GSE45291", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
petri_2 <- getGEO("GSE121239", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]


# Define expression & metadata ####
reduce_probes <- function(matr, genes){
  asplit <- split(1:nrow(matr), genes)
  expr <- do.call(rbind, lapply(asplit, function(x){
    subm <- matr[x, ]
    if(length(x) > 1){
      subm <- subm[sample(which(rowMeans(subm) == max(rowMeans(subm))), 1), ]
    }
    subm
  }))
  rownames(expr) <- names(asplit)
  
  expr
}

chiche_expr <- chiche@assayData$exprs
gene_info <- chiche@featureData@data
chiche_gene <- gene_info$`Gene symbol`
sample_info <- chiche@phenoData@data
treat <- sample_info$characteristics_ch1
treat <- gsub("group: ", "", treat, fixed = T)
chiche_treat <- gsub("Healthy control of SLE", "CTRL", treat, fixed = T)
chiche_expr <- reduce_probes(chiche_expr, chiche_gene)

banchereau_expr <- banchereau@assayData$exprs
gene_info <- banchereau@featureData@data
banchereau_gene <- gene_info$`Gene symbol`
sample_info <- banchereau@phenoData@data
treat <- sample_info$characteristics_ch1.10
treat <- gsub("disease state: ", "", treat, fixed = T)
banchereau_treat <- gsub("Healthy", "CTRL", treat, fixed = T)
banchereau_expr <- reduce_probes(banchereau_expr, banchereau_gene)

petri_1_expr <- petri_1@assayData$exprs
gene_info <- petri_1@featureData@data
petri_1_gene <- gene_info$`Gene symbol`
sample_info <- petri_1@phenoData@data
treat <- sample_info$`disease:ch1`
treat <- gsub("SLE (Systemic LUPUS Erythomatosus)", "SLE", treat, fixed = T)
treat <- gsub("Control", "CTRL", treat, fixed = T)
petri_1_treat <- treat
petri_1_expr <- reduce_probes(petri_1_expr, petri_1_gene)

petri_2_expr <- petri_2@assayData$exprs
gene_info <- petri_2@featureData@data
petri_2_gene <- gene_info$`Gene symbol`
sample_info <- petri_2@phenoData@data
treat <- sample_info$`disease state:ch1`
treat <- gsub("Systemic Lupus Erythematosus", "SLE", treat, fixed = T)
treat <- gsub("Healthy", "CTRL", treat, fixed = T)
petri_2_treat <- treat
petri_2_expr <- reduce_probes(petri_2_expr, petri_2_gene)


# Merge into 1 analysis ####
allgenes <- table(c(chiche_gene, banchereau_gene, petri_1_gene, petri_2_gene))
allgenes <- names(which(allgenes == 4))

merged_expr <- data.frame(
  chiche_expr[match(allgenes, rownames(chiche_expr)), ],
  banchereau_expr[match(allgenes, rownames(banchereau_expr)), ],
  petri_1_expr[match(allgenes, rownames(petri_1_expr)), ],
  petri_2_expr[match(allgenes, rownames(petri_2_expr)), ]
)

merged_treat <- c(
  chiche_treat,
  banchereau_treat,
  petri_1_treat,
  petri_2_treat
)

study <- c(
  rep("chiche", length(chiche_treat)),
  rep("banchereau", length(banchereau_treat)),
  rep("petri_1", length(petri_1_treat)),
  rep("petri_2", length(petri_2_treat))
)

other <- which(merged_treat %in% c(
  "Rheumatoid Arthiritis (DMARD-IR)", "Rheumatoid Arthiritis (TNF-IR)"))
merged_expr <- merged_expr[, -other]
merged_treat <- merged_treat[-other]
study <- study[-other]


# Run differential expression ####
design <- model.matrix(~ study  + merged_treat)
fit <- lmFit(merged_expr, design)
fit <- eBayes(fit)
top_table <- topTable(fit, coef = "merged_treatSLE", n = Inf)
top_table$gene <- rownames(top_table)


# Create volcano plot ####
ggplot(top_table, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05)) +
  scale_color_manual(values = c("grey", "darkred")) +
  labs(
    title = "SLE meta-analysis",
    y = "-log10(p-value)",
    x = "log2 fold change"
  ) +
  ggrepel::geom_text_repel(
    data = top_table[1:40, ],
    aes(label = gene),
  ) +
  theme_classic()

plot_gene <- function(gene){
  aframe <- data.frame(
    gene = as.numeric(merged_expr[gene, ]),
    study, merged_treat)
  
  ggplot(aframe, aes(merged_treat, gene, color = merged_treat)) +
    labs(
      title = "SLE meta-analysis",
      y = paste(gene, "expression levels"),
      x = "Condition", color = "Condition"
    ) +
    facet_wrap(~ study, nrow = 1, scales = "free_y") +
    geom_boxplot() +
    geom_jitter(width = 0.1) +
    scale_color_manual(values = c("grey", "darkred")) +
    theme_classic()
}
plot_gene("IFI27")
plot_gene("KLRB1")


# Save results in table ####
write.csv(
  top_table, file = "../outputs/merged_sle_diff_expr.csv",
  row.names = F
)


