library(Biobase); library(limma); library(ggplot2); library(ggrepel); library(pheatmap)

eset <- readRDS("data/GSE20692_mRNA_eset.rds")
tab  <- read.csv("results/tables/DEG_AgNP_vs_Control_batchAdj.csv")
titles <- pData(eset)$title
group  <- factor(ifelse(grepl("not exposed", titles), "Control",
                        ifelse(grepl("AgNPs", titles), "AgNP", "AgIon")),
                 levels = c("Control", "AgNP", "AgIon"))
rep_id <- factor(sub(".*rep([0-9]).*", "\\1", titles))

# --- Volcano ---
v <- tab
v$Status <- "Not significant"
v$Status[v$adj.P.Val < 0.05] <- "adj.P<0.05"
v$Status[v$adj.P.Val < 0.05 & abs(v$logFC) > 1] <- "adj.P<0.05 & |logFC|>1"
top <- head(v[order(v$adj.P.Val), ], 15)

p <- ggplot(v, aes(logFC, -log10(adj.P.Val), color = Status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top, aes(label = Symbol), size = 3, color = "black") +
  scale_color_manual(values = c("Not significant" = "grey70",
                                "adj.P<0.05" = "orange",
                                "adj.P<0.05 & |logFC|>1" = "firebrick")) +
  labs(title = "AgNP vs Control (Jurkat, 24 h)",
       x = "log2 fold change", y = "-log10(adjusted P-value)") +
  theme_bw()
ggsave("results/figures/06_volcano_AgNP_vs_Control.png", p, width = 8, height = 6, dpi = 300)

# --- Heatmap ۵۰ ژن برتر (اصلاح batch فقط برای نمایش) ---
top50 <- head(tab$ProbeID, 50)
mat <- removeBatchEffect(exprs(eset), batch = rep_id)[top50, ]
rownames(mat) <- make.unique(as.character(tab$Symbol[match(top50, tab$ProbeID)]))
annot <- data.frame(Group = group, row.names = colnames(mat))

png("results/figures/07_heatmap_top50.png", width = 900, height = 1100, res = 120)
pheatmap(mat, scale = "row", annotation_col = annot, fontsize_row = 7,
         main = "Top 50 genes: AgNP vs Control (batch-corrected for display)")
dev.off()