library(Biobase)
library(limma)

eset <- readRDS("data/GSE20692_mRNA_eset.rds")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# --- گروه‌ها ---
titles <- pData(eset)$title
group <- ifelse(grepl("not exposed", titles), "Control",
                ifelse(grepl("AgNPs", titles), "AgNP", "AgIon"))
group <- factor(group, levels = c("Control", "AgNP", "AgIon"))

# --- Design matrix (بدون intercept) ---
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

# --- برازش مدل ---
fit <- lmFit(exprs(eset), design)

# --- Contrastها ---
contrast_matrix <- makeContrasts(
  AgNP_vs_Control    = AgNP - Control,
  AgIon_vs_Control   = AgIon - Control,
  AgNP_vs_AgIon      = AgNP - AgIon,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# --- خلاصه‌ی تعداد ژن‌های معنی‌دار (adj.P < 0.05 و |log2FC| > 1) ---
res <- decideTests(fit2, method = "separate", adjust.method = "BH",
                   p.value = 0.05, lfc = 1)
summary(res)

# --- جدول نتایج هر contrast، همراه با نماد ژن ---
get_table <- function(coef_name) {
  tt <- topTable(fit2, coef = coef_name, number = Inf,
                 adjust.method = "BH", sort.by = "P")
  tt$ProbeID <- rownames(tt)
  tt$Symbol <- fData(eset)[tt$ProbeID, "Symbol"]
  tt
}

tab_np  <- get_table("AgNP_vs_Control")
tab_ion <- get_table("AgIon_vs_Control")
tab_dir <- get_table("AgNP_vs_AgIon")

head(tab_np[, c("Symbol", "logFC", "AveExpr", "P.Value", "adj.P.Val")], 15)

# --- اعتبارسنجی: ژن‌های متالوتیونین (مقاله‌ی اصلی گزارش کرده) ---
tab_np[grepl("^MT1|^MT2", tab_np$Symbol),
       c("Symbol", "logFC", "adj.P.Val")]

# --- ذخیره ---
write.csv(tab_np,  "results/tables/DEG_AgNP_vs_Control.csv",  row.names = FALSE)
write.csv(tab_ion, "results/tables/DEG_AgIon_vs_Control.csv", row.names = FALSE)
write.csv(tab_dir, "results/tables/DEG_AgNP_vs_AgIon.csv",    row.names = FALSE)