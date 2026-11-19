library(GEOquery)

gset <- getGEO("GSE20692", GSEMatrix = TRUE)

length(gset)                          # چند پلتفرم؟
gset[[1]]                             # خلاصه: تعداد نمونه و پلتفرم
table(pData(gset[[1]])$title)         # اسم نمونه‌ها
warnings()
options(timeout = 600)
options(download.file.method.GEOquery = "wininet")
library(GEOquery)
gset <- getGEO("GSE20692", GSEMatrix = TRUE)