library(Biobase)
library(ggplot2)
library(pheatmap)

eset <- readRDS("data/GSE20692_mRNA_eset.rds")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# --- تعریف گروه‌ها از روی عنوان نمونه‌ها ---
titles <- pData(eset)$title
group <- ifelse(grepl("not exposed", titles), "Control",
                ifelse(grepl("AgNPs", titles), "AgNP", "AgIon"))
group <- factor(group, levels = c("Control", "AgNP", "AgIon"))
table(group)           # باید 3 / 3 / 3 باشه
cbind(titles, as.character(group))   # چک کن گروه‌ها درست جا خوردن

ex <- exprs(eset)
cols <- c("grey60", "firebrick", "steelblue")[as.integer(group)]

# --- 1) Boxplot ---
png("results/figures/01_boxplot.png", width = 900, height = 600)
boxplot(ex, col = cols, las = 2, outline = FALSE,
        main = "Expression distribution per sample", ylab = "log2 intensity",
        names = paste0(group, 1:3))
dev.off()

# --- 2) Density plot ---
png("results/figures/02_density.png", width = 900, height = 600)
plot(density(ex[, 1]), col = cols[1], main = "Density of expression values",
     xlab = "log2 intensity", ylim = c(0, max(apply(ex, 2, function(x) max(density(x)$y)))))
for (i in 2:ncol(ex)) lines(density(ex[, i]), col = cols[i])
legend("topright", legend = levels(group),
       col = c("grey60", "firebrick", "steelblue"), lty = 1)
dev.off()

# --- 3) PCA ---
pca <- prcomp(t(ex), scale. = TRUE)
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     Group = group, Sample = colnames(ex))

p <- ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  labs(title = "PCA of samples",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_bw()
ggsave("results/figures/03_pca.png", p, width = 7, height = 5, dpi = 300)

# --- 4) Heatmap همبستگی نمونه‌ها ---
cor_mat <- cor(ex)
annot <- data.frame(Group = group, row.names = colnames(ex))
png("results/figures/04_sample_correlation.png", width = 800, height = 700)
pheatmap(cor_mat, annotation_col = annot, annotation_row = annot,
         main = "Sample-to-sample correlation")
dev.off()