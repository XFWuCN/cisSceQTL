#' Build gene-snp pairs metadata for eQTL-calling
#' The input snvmatrix from the user requires the encoding format of 0/1/2/3.
#' @param snvMatrix A genotype matrix where each row is one variant and each column is one sample, and the scoring method is 0/1/2/3.
#' @param genedata A gene expression matrix where each row is one gene and each column is one sample
#' @param scoring The user chooses whether to convert the counting method of the snvMatrix to -1/0/1, TRUE indicates conversion, and FALSE indicates no conversion, default is no conversion.
#' @param species The species that the user wants to select, human or mouse.
#' @param group
#' @param snv.number.of.cells
#' @param expression.min
#' @param expression.number.of.cells
#' @param gene_ids
#' @param cisDist
#' @importFrom org.Hs.eg.db
#' @importFrom biomaRt useMart getBM useEnsembl mapIds
#' @importFrom org.Mm.eg.db
#' @return
#' @export
#' @examples
#'
BuildMetadata <- function(snvMatrix,
                          genedata,
                          sparcity = FALSE,
                          species = "human",
                          group = "celltype",  # Group is celltype by default.
                          snv.number.of.cells = 30,
                          expression.min = 1,
                          expression.number.of.cells = 30,
                          gene_ids = NULL,
                          cisDist = NULL){

  if(class(genedata)[1] == "Seurat"){
    seurat_group <- genedata[[group]]
    unique_group <- unique(seurat_group[,1])
  }else if(is.matrix(genedata) | is.data.frame(genedata) | class(genedata)[1] == "dgCMatrix"){
    unique_group = "Matrix"
  }else{
    stop("请输入正确格式的数据！")
  }

  if(species == "human"){
    snv_dataset = "hsapiens_snp"
    gene_dataset = "hsapiens_gene_ensembl"
    OrgDb = org.Hs.eg.db
  }else if(species == "mouse"){
    snv_dataset = "mmusculus_snp"
    gene_dataset = "mmusculus_gene_ensembl"
    OrgDb = org.Mm.eg.db
  }else{
    stop("Please enter human or mouse.")
  }

  creat_snps_loc <- function(snv.list){

    # Get location for each SNV
    snp_mart <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = snv_dataset)

    snps_loc <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                      filters = "snp_filter",
                      values = snv.list,  # 指定染色体或snv_list
                      mart = snp_mart)
    # 把列名“chrom_start”改为“position”
    colnames(snps_loc)[which(colnames(snps_loc) == "chrom_start")] <- "position"
    rownames(snps_loc) <- snps_loc[, 1]
    return(snps_loc)
  }

  creat_gene_loc <- function(gene.list){

    # Get location for each gene
    gene_mart = useEnsembl(biomart="ensembl",
                           dataset="hsapiens_gene_ensembl")

    ensembls <- mapIds(OrgDb, keys = gene.list, keytype = "SYMBOL", column="ENSEMBL")
    ensembls = as.data.frame(ensembls)
    ensembls_id = unique(ensembls$ensembls)

    #ensembl_gene_id,external_gene_name
    gene_attributes=c('external_gene_name','chromosome_name','start_position','end_position')

    #获取gene的相关坐标信息
    gene_loc <- getBM(attributes = gene_attributes,
                      filters = "external_gene_name",
                      values = gene.list,  # 指定染色体或gene_list
                      mart = gene_mart)

    rownames(gene_loc) <- gene_loc[, 1]
    return(gene_loc)
  }



  check_snvList <- function(snv.list){
    if(grepl("^rs", snv.list[[1]][1])){  # chromosome start with "rs"
      creat_snps_loc(snv.list)
    } else if(grepl("\\d+:\\d+", snv.list[[1]][1])){  # chr:position
      # build an empty dataframe
      snps_df <- data.frame(refsnp_id = character(),
                            chr_name = character(),
                            position = numeric(),
                            stringsAsFactors = FALSE)
      snps_loc <- snps_df
      for (i in 1:length(snv.list)){
        snv_parts <- strsplit(snv.list[i], ":")[[1]]  # split chr_name and position
        snps_loc <- rbind(snps_loc, data.frame(refsnp_id = snv.list[i],
                                               chr_name = snv_parts[1],
                                               position = as.numeric(snv_parts[2])))
      }
      return(snps_loc)
    } else{
      stop("Error: SNV does not match expected format.")
    }
  }


  metadata <- data.frame()

  for(i in unique_group){
    print(i)
    if(i == "Matrix"){
      split_snvMatrix = snvMatrix
      split_expressionMatrix = genedata
    }else{
      logical_index <- genedata[[group]] == i
      split_seuratdata <- genedata[ ,logical_index]
      split_expressionMatrix <- as.matrix(split_seuratdata@assays$RNA@data)
      cell_list <- colnames(split_expressionMatrix)
      split_snvMatrix <- snvMatrix[, cell_list]
    }

    snv.list <- rownames(split_snvMatrix)
    gene.list <- rownames(split_expressionMatrix)
    cell.list <- colnames(split_expressionMatrix)

    # When users do not perform cis-eQTL filtering
    if(is.null(gene_ids) && is.null(cisDist)){
      matched_snps = snv.list
      matched_gene = gene.list
    }else if(!is.null(gene_ids) && is.null(cisDist)){
      matched_snps = snv.list
      # Check if gene_ids exist in gene.list
      not.in.list <- gene_ids[!gene_ids %in% gene.list]

      if (length(not.in.list) > 0) {
        # Prompt users with gene_ids that do not exist
        print("Invalid gene ID！")
      } else {
        matched_gene <- gene_ids
      }
    }else if(is.null(gene_ids) && !is.null(cisDist)){
      if (length(cisDist) != 2) {
        stop("cisDist should be specified as c(minDist, maxDist).")
      }
      if (cisDist[1] >= cisDist[2]) {
        stop("minDist should be smaller than maxDist.")
      }


      # 连接数据库，构建snps_loc和gene_loc

      # Get location for each SNV
      snps_loc <- check_snvList(snv.list)

      gene_loc <- creat_gene_loc(gene.list)


      # 遍历基因和SNV，匹配在基因范围内的SNV，输出gene-SNV pairs
      matched_gene <- c()
      matched_snps <- c()

      for (i in 1:nrow(gene_loc)) {
        # get range
        gene_start1 <- gene_loc$start_position[i] + cisDist[1]
        gene_end1 <- gene_loc$end_position[i] + cisDist[2]

        # 遍历每个SNV
        for (j in 1:nrow(snps_loc)) {
          # Retrieve the position information of the current SNV
          snp_chr <- snps_loc$chr_name[j]
          snp_pos <- snps_loc$position[j]

          # Determine if the current SNV is within the range of the gene
          if (snp_chr == gene_loc$chromosome_name[i] && snp_pos >= gene_start1 && snp_pos <= gene_end1) {
            # Add matching SNVs to the current gene's SNV list
            matched_snps <- c(matched_snps, snps_loc$refsnp_id[j])

            matched_gene <- c(matched_gene, gene_loc$ensembl_gene_id[i])

          }
        }

      }
      matched_snps <- unique(matched_snps)
      matched_gene <- unique(matched_gene)
    }else{ # 指定gene和cisDist
      if (length(cisDist) != 2) {
        stop("cisDist should be specified as c(minDist, maxDist).")
      }
      if (cisDist[1] >= cisDist[2]) {
        stop("minDist should be smaller than maxDist.")
      }

      # 连接数据库，构建snps_loc和gene_loc

      snps_loc <- check_snvList(snv.list)

      # 检查gene_ids是否存在于gene.list中
      not.in.list <- gene_ids[!gene_ids %in% gene.list]

      if (length(not.in.list) > 0) {
        # 提示用户不存在的gene_ids
        print("Invalid gene ID！")
      } else {
        gene.list <- gene_ids # 用户指定的genes
      }
      gene_loc <- creat_gene_loc(gene.list)


      # 遍历基因和SNV，匹配在基因范围内的SNV，输出gene-SNV pairs
      matched_gene <- c()
      matched_snps <- c()

      for (i in 1:nrow(gene_loc)) {
        # 获取范围
        gene_start1 <- gene_loc$start_position[i] + cisDist[1]
        gene_end1 <- gene_loc$end_position[i] + cisDist[2]

        # 遍历每个SNV
        for (j in 1:nrow(snps_loc)) {
          # 获取当前SNV的位置信息
          snp_chr <- snps_loc$chr_name[j]
          snp_pos <- snps_loc$position[j]

          # 判断当前SNV是否在基因的范围内
          if (snp_chr == gene_loc$chromosome_name[i] && snp_pos >= gene_start1 && snp_pos <= gene_end1) {
            # 将匹配的SNV添加到当前基因的SNV列表中
            matched_snps <- c(matched_snps, snps_loc$refsnp_id[j])

            matched_gene <- c(matched_gene, gene_loc$ensembl_gene_id[i])

          }
        }

      }
      matched_snps <- unique(matched_snps)
      matched_gene <- unique(matched_gene)
    }

    N <- length(matched_snps)

    if(sparcity == TRUE){

      # 构建空数据框，元素个数是SNV的个数N
      metadata_split <- data.frame(group = rep("", N),
                                   SNVid = rep("",N),
                                   Num_cells_ref = rep(0, N),
                                   Num_cells_alt = rep(0, N),
                                   Ref_cells = rep("", N),
                                   Alt_cells = rep("", N),
                                   GeneList = rep("", N),
                                   Num_gene = rep(0, N),
                                   CellList = rep("", N),
                                   stringsAsFactors = FALSE)

      useful_snv <- 0 # count snv number for final dataframe

      matched_snps[matched_snps == 0] <- -1
      matched_snps[matched_snps == 1] <- 0
      matched_snps[matched_snps == 2 | matched_snps == 3] <- 1


      for (snvid in matched_snps){
        cell.ref <- colnames(split_snvMatrix[snvid, split_snvMatrix[snvid,] == 0])
        cell.alt <- colnames(split_snvMatrix[snvid, split_snvMatrix[snvid,] == 1])
        if((length(cell.ref) > snv.number.of.cells) & (length(cell.alt) > snv.number.of.cells)){ # snv pass test
          # test gene further
          genelist <- c() # for saving genes for this snv
          for (gene in matched_gene){
            cell.ref.gene <- split_expressionMatrix[gene, cell.ref]
            cell.alt.gene <- split_expressionMatrix[gene, cell.alt]
            cell.ref.gene.valid.num = sum(cell.ref.gene > expression.min) # number of valid cells on this gene
            cell.alt.gene.valid.num = sum(cell.alt.gene > expression.min) # number of valid cells on this gene
            if ((cell.ref.gene.valid.num > expression.number.of.cells) & (cell.alt.gene.valid.num > expression.number.of.cells)){ # this gene pass the test
              genelist <- c(genelist, gene)
            }
          }
          if (length(genelist) > 0){ # have valid gene
            useful_snv = useful_snv + 1 # count valid snv number
            metadata_split[useful_snv, 'group'] <- i
            metadata_split[useful_snv, 'SNVid'] <- snvid
            metadata_split[useful_snv, 'Num_cells_ref'] <- length(cell.ref)
            metadata_split[useful_snv, 'Num_cells_alt'] <- length(cell.alt)
            metadata_split[useful_snv, 'Ref_cells'] <- paste(cell.ref, collapse = ',')
            metadata_split[useful_snv, 'Alt_cells'] <- paste(cell.alt, collapse = ',')
            metadata_split[useful_snv, 'GeneList'] <- paste(genelist, collapse = ',')
            metadata_split[useful_snv, 'Num_gene'] <- length(genelist)
            metadata_split[useful_snv, 'CellList'] <- paste(c(cell.ref,cell.alt), collapse = ',')
          }
        }
      }
      # remove nonsense rows
      metadata_split <- metadata_split[1:useful_snv, ]
      metadata <- rbind(metadata, metadata_split)
    }else if(sparcity == FALSE){

      # 构建空数据框，元素个数是SNV的个数N
      metadata_split <- data.frame(group = rep("", N),
                                   SNVid = rep("", N),
                                   Num_cells_AA = rep(0, N),
                                   Num_cells_Aa = rep(0, N),
                                   Num_cells_aa = rep(0, N),
                                   AA_cells = rep("", N),
                                   Aa_cells = rep("", N),
                                   aa_cells = rep("", N),
                                   GeneList = rep("", N),
                                   Num_gene = rep(0, N),
                                   CellList = rep("", N),
                                   stringsAsFactors = FALSE)

      useful_snv <- 0 # count snv number for final dataframe

      for (snvid in matched_snps){
        cell.AA <- colnames(split_snvMatrix[, split_snvMatrix[snvid,] == 1])
        cell.Aa <- colnames(split_snvMatrix[, split_snvMatrix[snvid,] == 3])
        cell.aa <- colnames(split_snvMatrix[, split_snvMatrix[snvid,] == 2])
        if((length(cell.AA) > snv.number.of.cells) & (length(cell.Aa) > snv.number.of.cells) & (length(cell.aa) > snv.number.of.cells)){ # snv pass test
          # test gene further
          genelist <- c() # for saving genes for this snv
          for (gene in matched_gene){
            cell.AA.gene <- split_expressionMatrix[gene, cell.AA]
            cell.Aa.gene <- split_expressionMatrix[gene, cell.Aa]
            cell.aa.gene <- split_expressionMatrix[gene, cell.aa]

            cell.AA.gene.valid.num = sum(cell.AA.gene > expression.min)
            cell.Aa.gene.valid.num = sum(cell.Aa.gene > expression.min)
            cell.aa.gene.valid.num = sum(cell.aa.gene > expression.min)
            if ((cell.AA.gene.valid.num > expression.number.of.cells) & (cell.Aa.gene.valid.num > expression.number.of.cells) & (cell.aa.gene.valid.num > expression.number.of.cells)){ # this gene pass the test
              genelist <- c(genelist, gene)
            }
          }
          if (length(genelist) > 0){ # have valid gene
            useful_snv = useful_snv + 1 # count valid snv number
            #metadata_split[useful_snv,'rsID'] <- snv_id
            metadata_split[useful_snv, 'group'] <- i
            metadata_split[useful_snv, 'SNVid'] <- snvid
            metadata_split[useful_snv, 'Num_cells_AA'] <- length(cell.AA)
            metadata_split[useful_snv, 'Num_cells_Aa'] <- length(cell.Aa)
            metadata_split[useful_snv, 'Num_cells_aa'] <- length(cell.aa)
            metadata_split[useful_snv, 'AA_cells'] <- paste(cell.AA, collapse = ',')
            metadata_split[useful_snv, 'Aa_cells'] <- paste(cell.Aa, collapse = ',')
            metadata_split[useful_snv, 'aa_cells'] <- paste(cell.aa, collapse = ',')
            metadata_split[useful_snv, 'GeneList'] <- paste(genelist, collapse = ',')
            metadata_split[useful_snv, 'Num_gene'] <- length(genelist)
            metadata_split[useful_snv, 'CellList'] <- paste(c(cell.AA,cell.Aa,cell.aa), collapse = ',')
          }
        }
      }
      # remove nonsense rows
      metadata_split <- metadata_split[1:useful_snv, ]
      metadata <- rbind(metadata, metadata_split)
    }else{
      stop("sparcity can only be selected as 'TRUE' or 'FALSE'")
    }
  }
  return(metadata)
}
