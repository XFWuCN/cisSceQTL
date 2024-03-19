#' Title
#'
#' @param expressionMatrix 基因表达矩阵
#' @param method 进行标准化的方法，此处有四个方法可选
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @importFrom limma normalizeBetweenArrays
#' @importFrom sctransform vst
#' @return 返回标准化后的数据
#' @export
#'
#' @examples
#'
#'
#'

Normalize <- function(expressionMatrix, method = "sctransform") {
  if (!method %in% c("CPM", "TPM", "DESeq", "limma", "sctransform")) {
    stop("Invalid method. Please choose from 'CPM', 'TPM', 'DESeq','limma' or 'sctransform'.")
  }

  # 计算每行的总计数，即每个基因在所有细胞的表达总量
  rowsum = apply(expressionMatrix, 1, sum)
  # 将没有在这个矩阵的任何细胞中表达的基因移除
  expressionMatrix = expressionMatrix[rowsum!=0,]

  if (method == "CPM") {
    # CPM normalization for non-zero values
    total_counts <- colSums(expressionMatrix)
    normalized_data <- sweep(expressionMatrix, 2, total_counts, "/") * 1e6
  } else if (method == "TPM") {
    # TPM normalization for non-zero values
    library_size <- colSums(expressionMatrix)
    normalized_data <- expressionMatrix / library_size * 1e6
  } else if (method == "DESeq") {
    # DESeq normalization
    # 先将表达矩阵转换为DESeqDataSet对象
    # DESeq normalization
    library("DESeq2")
    sample.df <- colnames(expressionMatrix)
    sample.df <- as.data.frame(sample.df)  # 这一步使得sample.df有维度，否则下一步行不通
    dds <- DESeqDataSetFromMatrix(countData = expressionMatrix,
                                  colData = sample.df,
                                  design = ~1)
    dds <- DESeq(dds)
    normalized_data <- counts(dds, normalized = TRUE)
  } else if (method == "limma") {
    # limma normalization
    library("limma")
    normalized_data <- normalizeBetweenArrays(expressionMatrix, method="quantile")
  } else if (method == "sctransform") {
    # sctransform normalization
    normalized_data <- as.matrix(expressionMatrix)
    normalized_data <- sctransform::vst(normalized_data)$y
  }

  # 输出标准化结果信息
  cat(paste("Normalization completed using method:", method, "\n"))
  cat("Dimensions of normalized data:", dim(normalized_data), "\n")
  cat("Summary statistics of samples:\n")
  colSum <- colSums(normalized_data)
  colMean <- colMeans(normalized_data)

  for (i in seq_along(colSum)) {
    cat(paste("\nSample", i, ":\n"," Total count =", colSum[i], "\n","Mean =",colMean[i],"\n"))
  }

  return(normalized_data)
}

