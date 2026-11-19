library(GEOquery)

# ساخت پوشه‌های لازم (اگه وجود ندارن)
dir.create("data", showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# دانلود داده از GEO
options(timeout = 600)
gset <- getGEO("GSE20692", GSEMatrix = TRUE)

# نگه داشتن فقط mRNA (پلتفرم Illumina GPL6883)
eset <- gset[["GSE20692-GPL6883_series_matrix.txt.gz"]]

# ذخیره برای اسکریپت‌های بعدی
saveRDS(eset, "data/GSE20692_mRNA_eset.rds")