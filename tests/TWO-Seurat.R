
# snvMatrix = REF/ALT, genedata = seuratdata ------------------------------
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


# build metadata ----------------------------------------------------------

# 二分类的metadata
m2 <- BuildMetadata(snvMatrix = snvMatrix2,
                    genedata = seuratdata2,
                    sparcity = TRUE,
                    species = "human",
                    group = "celltype",
                    snv.number.of.cells = 0,
                    expression.min = 0,
                    expression.number.of.cells = 0,
                    gene_ids = NULL,
                    cisDist = NULL)



# call sceQTL -------------------------------------------------------------

# eqtl_zinb2出现MLE收敛警告，有70个snv, 没有GMP的细胞类型，全是CMP

eqtl_possion2 <- eQTLsc(seuratdata2, metadata = m2, sparcity = T,
                        p.adjust.method = "BH", useModel = "possion", p.adjust.Threshold = 1)


# choose a snp-gene pair to visualize -------------------------------------

p1 <- visualize.sample(eqtl = eqtl_possion2, "CMP", "1:632647", "CNN2", seuratdata2, m2, sparcity = TRUE, plottype='boxplot', removeoutlier = TRUE)
p2 <- visualize.sample(eqtl = eqtl_possion2, "GMP", "1:632647", "CNN2", seuratdata2, m2, sparcity = TRUE, plottype='boxplot', removeoutlier = TRUE)

grid.arrange(p1, p2, ncol = 2)
# 添加标题
grid.text("boxplot of CNN2 and 1:632647", x = unit(0.5, "npc"), y = unit(0.95, "npc"), just = "center", gp = gpar(fontsize = 26))
