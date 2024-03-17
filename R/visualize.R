#' Title
#'
#' @param eqtl
#' @param group
#' @param SNVid
#' @param Geneid
#' @param genedata
#' @param metadata
#' @param sparcity
#' @param plottype
#' @param removeoutlier
#'
#' @return
#' @export
#'
#' @examples
visualize.sample <- function(eqtl, group = "celltype", SNVid, Geneid, genedata, metadata, sparcity = FALSE, plottype='boxplot', removeoutlier = TRUE){

  if(class(genedata)[1] == "Seurat"){
    expressionMatrix = as.data.frame(genedata@assays$RNA@data)
    # 基因表达矩阵的处理可以加在这里
  }else if(is.matrix(genedata) | is.data.frame(genedata) | class(genedata)[1] == "dgCMatrix"){
    expressionMatrix = genedata
  }else{
    stop("Please enter the data in the correct format！")
  }

  result <- eqtl[(eqtl$group == group)&(eqtl$SNVid == SNVid)&(eqtl$Geneid == Geneid),]




  if(sparcity == TRUE){
    #collecting data from metadata
    SNVid <- result$SNVid
    geneid <- result$Geneid
    group <- result$group
    adjust.pvalue <- result$adjusted_pvalue
    ref_cells <- metadata[metadata$SNVid==SNVid, 'Ref_cells']
    ref_cells <- stringr::str_split(ref_cells, ",")[[1]]
    alt_cells <- metadata[metadata$SNVid==SNVid, 'Alt_cells']
    alt_cells <- stringr::str_split(alt_cells, ",")[[1]]


    if (removeoutlier) {
      sample_gene = expressionMatrix[geneid, ]
      non_zero_mask = sample_gene != 0
      sample_no_zero = sample_gene[, colSums(non_zero_mask) > 0]
      med = median(unlist(sample_no_zero))
      mad = mad(unlist(sample_no_zero))
      filter = sample_no_zero[, colMeans(sample_no_zero) < med + 4 * mad]
      ref_cells = intersect(ref_cells, colnames(filter))
      alt_cells = intersect(alt_cells, colnames(filter))
    }

    counts_Ref <- unlist(expressionMatrix[geneid, ref_cells])
    counts_Alt <- unlist(expressionMatrix[geneid, alt_cells])

    # building data frame
    df <- data.frame(expression = c(counts_Ref,counts_Alt),
                     snp = c(rep("REF",length(counts_Ref)), rep("ALT",length(counts_Alt))))
    df$snp <- factor(df$snp, levels = c("REF", "ALT"))
    # title <- paste(plottype, "of", geneid , "and", SNVid, "in", group )
    title <- paste(group)


  }else if(sparcity == FALSE){
    SNVid <- result$SNVid
    geneid <- result$Geneid
    group <- result$group
    adjust.pvalue <- result$adjusted_pvalue
    AA_cells <- metadata[metadata$SNVid==SNVid, 'AA_cells']
    AA_cells <- stringr::str_split(AA_cells, ",")[[1]]
    Aa_cells <- metadata[metadata$SNVid==SNVid, 'Aa_cells']
    Aa_cells <- stringr::str_split(Aa_cells, ",")[[1]]
    aa_cells <- metadata[metadata$SNVid==SNVid, 'aa_cells']
    aa_cells <- stringr::str_split(aa_cells, ",")[[1]]

    if (removeoutlier){
      sample_gene = expressionMatrix[geneid,]
      sample_no_zero = sample_gene[,sample_gene != 0]
      med = median(unlist(sample_no_zero))
      mad = mad(unlist(sample_no_zero))
      filter = sample_no_zero[,sample_no_zero < med+4*mad]
      AA_cells = intersect(AA_cells, colnames(filter))
      Aa_cells = intersect(Aa_cells, colnames(filter))
      aa_cells = intersect(aa_cells, colnames(filter))
    }

    counts_AA <- unlist(expressionMatrix[geneid, AA_cells])
    counts_Aa <- unlist(expressionMatrix[geneid, Aa_cells])
    counts_aa <- unlist(expressionMatrix[geneid, aa_cells])

    # building data frame
    df <- data.frame(expression = c(counts_AA,counts_Aa,counts_aa),
                     snp = c(rep("AA",length(counts_AA)), rep("Aa",length(counts_Aa)), rep("aa",length(counts_aa))))
    df$snp <- factor(df$snp, levels = c("AA", "Aa", "aa"))
    # title <- paste(plottype, "of", geneid , "and", SNVid, "in", group)
    title <- paste(group)
  }else{
    stop("sparcity can only be selected as 'TRUE' or 'FALSE'")
  }


  drawboxplot <- function(df){
    ggplot(df, aes(x = factor(df$snp),
                   y = df$expression,
                   fill = factor(df$snp)))+
      geom_boxplot(alpha=0.3)+
      theme_ipsum()+
      ggtitle(title)+
      scale_fill_brewer(palette="Dark2")+
      labs(x = "Group", y = "Expression")+
      theme(plot.title = element_text(size = 14),
            axis.line = element_line(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))
  }

  drawviolinplot <- function(df){
    ggplot(df, aes(x = factor(df$snp),
                   y = df$expression,
                   fill = factor(df$snp)))+
      geom_violin(alpha=0.5)+
      theme_ipsum()+
      geom_boxplot(width = 0.12, outlier.shape = NA, position = position_dodge(0.55)) +  # 显示箱线图
      ggtitle(title)+
      scale_fill_brewer(palette="Dark2")+
      labs(x = "Group", y = "Expression")+
      theme(plot.title = element_text(size = 14),
            axis.line = element_line(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))

  }

  drawhistplot <- function(df){
    ggplot(df, aes(df$expression,fill = factor(df$snp)))+
      geom_histogram(binwidth=40,color="white")+
      scale_fill_brewer(palette = "Pastel1")+
      facet_grid(df$snp ~ .,scales="free_y")+
      labs(title = "Expression Distribution by SNP Group",
           x = "Gene Expression",
           y = "Frequency")+
      theme_minimal()+
      ggtitle(title)+
      theme(legend.position = "top")+
      guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))
  }

  if(plottype=='violin'){
    drawviolinplot(df)
  }else if(plottype=='boxplot'){
    drawboxplot(df)
  }else if(plottype=='histplot'){
    drawhistplot(df)
  }


}
