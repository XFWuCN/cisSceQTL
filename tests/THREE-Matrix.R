
# snvMatrix = AA/aa/Aa, genedata = Matrix ------------------------------
library(digest)
library(cli)
library(GPLVM)
library(gaussianProcess)
library(nloptwrap)
library(mixtools)
library(Matrix)
library(lme4)
library(stats4)
library(glmmTMB)
library(stringr)
library(MASS)
library(bbmle)
library(splines)
library(gamlss.data)
library(gamlss.dist)
library(nlme)
library(parallel)
library(gamlss)
library(stats)
library(miscTools)
library(maxLik)
library(Matrix)
library(VGAM)
library(pscl)

library(org.Hs.eg.db)

library(ggplot2)
library(hrbrthemes)
library(tidyverse)
library(ggthemes)
library(showtext)
library(gridExtra)
library(grid)

# install.packages("devtools")
library(devtools)
install_github("XFWuCN/cisSceQTL")


# library test data -------------------------------------------------------

seuratdata <- readRDS("data/test_seurat_data.rds")
snvMatrix <- read.table("data/test_snp_data.txt", stringsAsFactors = F)
column_names <- colnames(snvMatrix)
new_column_names <- sapply(strsplit(column_names, ".1"), function(x) x[1])
intersection_string <- paste0(new_column_names, "-1")
colnames(snvMatrix) <- intersection_string

cell_list <- colnames(seuratdata)
random_cells <- sample(cell_list, 100)
seuratdata2 <- seuratdata[, random_cells]
snvMatrix2 <- snvMatrix[,random_cells]
gene_list <- rownames(seuratdata)
random_genes <- sample(gene_list, 20)
# seurat
seuratdata2 <- seuratdata2[random_genes, ]
# 处理seurat数据里的小数点
seuratdata2@assays$RNA@data <- floor(seuratdata2@assays$RNA@data * 1000)

snvMatrix2 <- snvMatrix2[1:100, ]

# Matrix
expressionMatrix2 = as.data.frame(seuratdata2@assays$RNA@data)


# build metadata ----------------------------------------------------------

# 三分类的的仅提供基因表达矩阵的metadata
mm3 <- BuildMetadata(snvMatrix = snvMatrix2,
                     genedata = expressionMatrix2,
                     sparcity = FALSE,
                     species = "human",
                     snv.number.of.cells = 0,
                     expression.min = 0,
                     expression.number.of.cells = 0,
                     gene_ids = NULL,
                     cisDist = NULL)



# call sceQTL -------------------------------------------------------------

# 仅提供基因表达矩阵
e_possion3 <- eQTLsc(expressionMatrix2, metadata = mm3, sparcity = F,
                     p.adjust.method = "BH", useModel = "possion", p.adjust.Threshold = 1)

e_linear3 <- eQTLsc(expressionMatrix2, metadata = mm3, sparcity = F,
                     p.adjust.method = "BH", useModel = "linear", p.adjust.Threshold = 1)

e_zinb3 <- eQTLsc(expressionMatrix2, metadata = mm3, sparcity = F,
                     p.adjust.method = "BH", useModel = "zinb", p.adjust.Threshold = 1)




# choose a snp-gene pair to visualize -------------------------------------

visualize.sample(eqtl = e_possion3, "Matrix", SNVid = "1:632647", Geneid = "CNN2", genedata = expressionMatrix2, metadata = mm3, sparcity = FALSE, plottype='boxplot', removeoutlier = TRUE)
