# نصب BiocManager (مدیر پکیج‌های Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# نصب پکیج‌های اصلی از Bioconductor
BiocManager::install(c("GEOquery", "limma", "Biobase"))

# نصب پکیج‌های رسم نمودار و کار با داده از CRAN
install.packages(c("ggplot2", "pheatmap", "dplyr"))

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
install.packages(c("ggplot2", "pheatmap", "dplyr"))
install.packages("dplyr")
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
