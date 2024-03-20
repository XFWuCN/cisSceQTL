#' eQTLsc: Uncover single-cell eQTLs exclusively using scRNA-seq data
#'
#'
#' A function designed to identify eQTLs from scRNA-seq data,
#' utilizing the eQTLsingle_build_metadata function to generate SNP-gene pair information.
#' 这个函数用于从单细胞RNA测序数据中识别eQTLs，并利用BuildMetadata函数生成SNP-gene对信息。
#' @param expressionMatrix 基因表达矩阵
#' @param metadata gene和snp的pairs
#' @param dist 指定使用负二项分布作为数据的分布类型
#' @param EM 指定使用EM算法进行参数估计
#' @param p.adjust.method p值校验
#' @param theta 表示模型中过度离散的参数，用于描述计数数据的离散程度
#' @param mu 表示负二项分布的均值，用于描述非零计数的平均水平
#' @param size 表示负二项分布的大小参数，用于描述数据的离散程度
#' @param prob 表示零膨胀分布的概率参数，用于描述计数为零的概率
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb fitdistr
#' @importFrom VGAM dzinegbin
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom pscl zeroinfl
#' @importFrom stats p.adjust pchisq plogis
#' @importFrom glmmTMB glmmTMB
#' @return A dataframe, each row describes eQTL discovering result of a snv-gene pair
#' 返回一个数据框，每一行描述了一个SNP-gene对的eQTL发现结果
#' @examples
#' @export


eQTLsc <- function(genedata, snvMatrix = NULL, metadata, sparcity = FALSE, p.adjust.method = "bonferroni", useModel = "ZINB", p.adjust.Threshold = 0.05){

  if(class(genedata)[1] == "Seurat"){
    expressionMatrix = as.data.frame(genedata@assays$RNA@data)
    # 基因表达矩阵的处理可以加在这里
  }else if(is.matrix(genedata) | is.data.frame(genedata) | class(genedata)[1] == "dgCMatrix"){
    expressionMatrix = genedata
  }else{
    stop("请输入正确格式的数据！")
  }


  # 格式校正 -----------------------------------------------
  if (sum(is.na(expressionMatrix)) > 0){
    stop("NA detected in 'expressionMatrix'");gc();
  }
  if (sum(expressionMatrix < 0) > 0) {
    stop("Negative value detected in 'expressionMatrix'")
  }
  if (all(expressionMatrix == 0)) {
    stop("All elements of 'expressionMatrix' are zero")
  }

  LINEAR <- function(expressionMatrix, snvMatrix, metadata, sparcity = FALSE, p.adjust.method = "bonferroni", p.adjust.Threshold = 0.05){

    result_all <- data.frame()

    unique_group <- unique(metadata$group)

    for(k in unique_group){

      result <- data.frame(
        SNVid = character(),
        group = character(),
        Geneid = character(),
        pvalue = double(),
        adjusted_pvalue = double(),
        Remark = character(),
        stringsAsFactors=FALSE)

      metadata_split1 <- metadata$group == k
      metadata_split <- metadata[metadata_split1, ]


      if(sparcity == FALSE){
        for(i in 1:dim(metadata_split)[1]){
          snvid <- metadata_split[i, "SNVid"]
          snv_mat <- snvMatrix[snvid, ]
          snv_mat <- t(snv_mat)

          genes <- metadata_split[i, 'GeneList']
          genes <- stringr::str_split(genes, ",")[[1]]
          for(j in 1:length(genes)){
            gene_id <- genes[j]
            gene_mat <- expressionMatrix[gene_id, ]
            gene_mat <- t(gene_mat)

            lmodel = lm( gene_mat ~ snv_mat);

            lmout = summary(lmodel)$coefficients[2, "Pr(>|t|)"]
            new_row <- data.frame(SNVid = snvid, group = k, Geneid = genes[j], pvalue = lmout)
            result <- rbind(result, new_row)
          }

        }
      }else if(sparcity == TRUE){

        snvMatrix[snvMatrix == 0] <- -1
        snvMatrix[snvMatrix == 1] <- 0
        snvMatrix[snvMatrix == 2 | snvMatrix == 3] <- 1

        for(i in 1:dim(metadata_split)[1]){
          snvid <- metadata_split[i, "SNVid"]
          snv_mat <- snvMatrix[snvid, ]
          snv_mat <- t(snv_mat)

          genes <- metadata_split[i, 'GeneList']
          genes <- stringr::str_split(genes, ",")[[1]]
          for(j in 1:length(genes)){
            gene_id <- genes[j]
            gene_mat <- expressionMatrix[gene_id, ]
            gene_mat <- t(gene_mat)

            lmodel = lm( gene_mat ~ snv_mat);

            lmout = summary(lmodel)$coefficients[2, "Pr(>|t|)"]
            new_row <- data.frame(SNVid = snvid, group = k, Geneid = genes[j], pvalue = lmout)
            result <- rbind(result, new_row)
          }

        }
      }else{
        stop("sparcity can only be selected as 'TRUE' or 'FALSE'")
      }


      # 加p值矫正方法
      if (!p.adjust.method %in% c("bonferroni", "holm", "hochberg", "hommel", "BH")) {
        stop("Invalid p-adjusted method.
           Please choose from 'bonferroni', 'holm', 'hochberg', 'hommel', or'fdr or BH'.")
      }

      # adjust p-value
      result[,"adjusted_pvalue"] <- p.adjust(result[,"pvalue"], method = "BH")
      # order adjust p-value
      result <- result[order(result[,"adjusted_pvalue"]),]
      # 重新编号
      rownames(result) <- NULL
      # 根据阈值过滤SNV-gene pairs
      result <- result[result$adjusted_pvalue <= p.adjust.Threshold, ]
      result_all <- rbind(result_all, result)

    }
    return(result_all)
  }




  ZINB <- function(expressionMatrix, snvMatrix = NULL, metadata, sparcity = FALSE, p.adjust.method = "bonferroni", p.adjust.Threshold = 1e-5){
    # 这里的参数千万不能回车！

    # p-value correction methods
    if (!p.adjust.method %in% c("bonferroni", "holm", "hochberg", "hommel", "BH")) {
      stop("Invalid p-adjusted method.
           Please choose from 'bonferroni', 'holm', 'hochberg', 'hommel', or'fdr or BH'.")
    }

    unique_group <- unique(metadata$group)

    result_all <- data.frame()

    for(j in unique_group){
      if(sparcity == TRUE){
        eQTLcalling <- function(i){
          if (i %% 100 == 0) {
            gc()
          }

          # get data from metadata : rsID,snvid,cells,ref_cells,alt_cells,genes ------
          snvid <- metadata_split[i, "SNVid"]
          cells_character <- metadata_split[i, 'CellList']
          cells <- stringr::str_split(cells_character, ",")[[1]]
          ref_cells_character <- metadata_split[i, 'Ref_cells']
          ref_cells <- stringr::str_split(ref_cells_character, ",")[[1]]
          alt_cells_character <- metadata_split[i, 'Alt_cells']
          alt_cells <- stringr::str_split(alt_cells_character, ",")[[1]]
          genes <- metadata_split[i, 'GeneList']
          genes <- stringr::str_split(genes, ",")[[1]]
          gene.cnt <- 0

          # snp ----------------------------------------------------------------------
          results_snp <- data.frame(
            SNVid = character(),
            group = character(),
            Geneid = character(),
            sample_size_1 = integer(),
            sample_size_2 = integer(),
            theta_1 = double(),
            theta_2 = double(),
            mu_1 = double(),
            mu_2 = double(),
            size_1 = double(),
            size_2 = double(),
            prob_1 = double(),
            prob_2 = double(),
            total_mean_1 = double(),
            total_mean_2 = double(),
            foldChange = double(),
            chi = double(),
            pvalue = double(),
            adjusted_pvalue = double(),
            Remark = character(),
            stringsAsFactors = FALSE)

          # gene ---------------------------------------------------------------------
          for (gene in genes){
            gene.cnt <- gene.cnt + 1
            # gene expression for ref group
            counts_1 <- unlist(expressionMatrix[gene, ref_cells])  # 从表达矩阵中选择指定基因在ref_cells中的counts值
            # gene expression for alt group
            counts_2 <- unlist(expressionMatrix[gene, alt_cells])
            results_gene <- data.frame(group = j,
                                       SNVid = snvid,
                                       Geneid = gene,
                                       sample_size_1 = length(counts_1),
                                       sample_size_2 = length(counts_2),
                                       theta_1 = NA,
                                       theta_2 = NA,
                                       mu_1 = NA,
                                       mu_2 = NA,
                                       size_1 = NA,
                                       size_2 = NA,
                                       prob_1 = NA,
                                       prob_2 = NA,
                                       total_mean_1 = NA,
                                       total_mean_2 = NA,
                                       foldChange = NA,
                                       chi = NA,
                                       pvalue = NA,
                                       adjusted_pvalue = NA,
                                       Remark = NA,
                                       stringsAsFactors = FALSE)

            # calculate fold change --------------------------------------------------
            totalMean_1 <- mean(counts_1)
            totalMean_2 <- mean(counts_2)
            foldChange  <- totalMean_1/totalMean_2

            # filling data -----------------------------------------------------------
            results_gene[1,"total_mean_1"] <- totalMean_1
            results_gene[1,"total_mean_2"] <- totalMean_2
            results_gene[1,"foldChange"]   <- foldChange


            # build  model -------------------------------------------------------

            build_model <- function(counts) {

              # 用gamlssML拟合ZINB模型
              options(show.error.messages = FALSE)  # 在代码运行时不显示错误消息
              zinb_gamlssML <- try(gamlssML(counts, family = "ZINBI"), silent = TRUE)
              options(show.error.messages = TRUE)  # 在运行代码时显示错误消息

              if('try-error' %in% class(zinb_gamlssML)){  # 如果建模失败
                print("MLE of ZINB failed! Please choose another model")
                return(list(theta = NA, mu = NA, size = NA, prob = NA))
              }else{
                # 拟合成功，获取参数
                zinb <- zinb_gamlssML
                theta <- zinb$nu
                mu <- zinb$mu
                size <- 1 / zinb$sigma
                prob <- size / (size + mu)
              }
              return(list(theta = theta, mu = mu, size = size, prob = prob))
            }

            build_model(counts_1) # 调用qtl-calling函数，传参
            build_model(counts_2)


            # estimate Ref params ----------------------------------------------------
            params_1 <- build_model(counts_1)
            theta_1 <- params_1[['theta']]
            mu_1 <- params_1[['mu']]
            size_1 <- params_1[['size']]
            prob_1 <- params_1[['prob']]
            # print(paste("counts_1",counts_1,"theta_1:", theta_1, "mu_1:", mu_1, "size_1:", size_1, "prob_1:", prob_1))


            #estimate Alt params -----------------------------------------------------
            params_2 <- build_model(counts_2)
            theta_2 <- params_2[['theta']]
            mu_2 <- params_2[['mu']]
            size_2 <- params_2[['size']]
            prob_2 <- params_2[['prob']]
            #print(paste("counts_2",counts_2,"theta_2:", theta_2, "mu_2:", mu_2, "size_2:", size_2, "prob_2:", prob_2))

            #combine Ref and Alt -----------------------------------------------------
            params_combined <- build_model(c(counts_1, counts_2))
            theta_res <- params_combined[['theta']]
            mu_res <- params_combined[['mu']]
            size_res <- params_combined[['size']]
            prob_res <- params_combined[['prob']]
            #print(paste("theta_res:", theta_res, "mu_res:", mu_res, "size_res:", size_res, "prob_res:", prob_res))


            # calculate p-value ---------------------------------------
            # 定义函数logL,计算给定参数下的对数似然函数值
            logL <- function(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2){
              # log-likelihood for count1 under parameter
              logL_1 <- sum(dzinegbin(counts_1, size = size_1, prob = prob_1, pstr0 = theta_1, log = TRUE))
              # log-likelihood for count2 under parameter
              logL_2 <- sum(dzinegbin(counts_2, size = size_2, prob = prob_2, pstr0 = theta_2, log = TRUE))
              logL <- logL_1 + logL_2
              logL
            }
            # 用上面定义的函数计算两个模型的对数似然函数值
            logL_1 <- logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2)
            logL_2 <- logL(counts_1, theta_res, size_res, prob_res, counts_2, theta_res, size_res, prob_res)
            chi <- logL_1 - logL_2
            pvalue <- 1 - pchisq(2 * chi , df = 3)


            # filling data into data frame ------------------------------

            results_gene[1,"theta_1"] <- theta_1
            results_gene[1,"theta_2"] <- theta_2
            results_gene[1,"mu_1"] <- mu_1
            results_gene[1,"mu_2"] <- mu_2
            results_gene[1,"size_1"] <- size_1
            results_gene[1,"size_2"] <- size_2
            results_gene[1,"prob_1"] <- prob_1
            results_gene[1,"prob_2"] <- prob_2
            results_gene[1,"chi"] <- chi
            results_gene[1,"pvalue"] <- pvalue

            results_snp <- rbind(results_snp, results_gene)
          }

          # return
          return(results_snp)

          # eQTLcalling end here
        }


        # final result
        result <- data.frame(
          SNVid = character(),
          group = character(),
          Geneid = character(),
          sample_size_1 = integer(),
          sample_size_2 = integer(),
          theta_1 = double(),
          theta_2 = double(),
          mu_1 = double(),
          mu_2 = double(),
          size_1 = double(),
          size_2 = double(),
          prob_1 = double(),
          prob_2 = double(),
          total_mean_1 = double(),
          total_mean_2 = double(),
          foldChange = double(),
          chi = double(),
          pvalue = double(),
          adjusted_pvalue = double(),
          Remark = character(),
          stringsAsFactors=FALSE)

        metadata_split1 <- metadata$group == j
        metadata_split <- metadata[metadata_split1, ]

        #process visible
        message("start calling eQTL")
        for(i in 1:dim(metadata_split)[1]){
          message(paste0("processing ", (i / dim(metadata_split)[1])*100,"%"))
          # 调用eQTLcalling函数传参给eQTLsc function
          result <- rbind(result, eQTLcalling(i))
        }
        message("finished!")


        #change column names
        colnames(result)[colnames(result) == 'sample_size_1'] <- 'sample_size_Ref'
        colnames(result)[colnames(result) == 'sample_size_2'] <- 'sample_size_Alt'
        colnames(result)[colnames(result) == 'theta_1'] <- 'theta_Ref'
        colnames(result)[colnames(result) == 'theta_2'] <- 'theta_Alt'
        colnames(result)[colnames(result) == 'mu_1'] <- 'mu_Ref'
        colnames(result)[colnames(result) == 'mu_2'] <- 'mu_Alt'
        colnames(result)[colnames(result) == 'size_1'] <- 'size_Ref'
        colnames(result)[colnames(result) == 'size_2'] <- 'size_Alt'
        colnames(result)[colnames(result) == 'prob_1'] <- 'prob_Ref'
        colnames(result)[colnames(result) == 'prob_2'] <- 'prob_Alt'
        colnames(result)[colnames(result) == 'total_mean_1'] <- 'total_mean_Ref'
        colnames(result)[colnames(result) == 'total_mean_2'] <- 'total_mean_Alt'
      }else if(sparcity == FALSE){
        eQTLcalling <- function(i){
          if (i %% 100 == 0) {
            gc()
          }

          snvid <- metadata_split[i, "SNVid"]
          cells_character <- metadata_split[i, 'CellList']
          cells <- stringr::str_split(cells_character, ",")[[1]]
          AA_cells_character <- metadata_split[i, 'AA_cells']
          AA_cells <- stringr::str_split(AA_cells_character, ",")[[1]]
          Aa_cells_character <- metadata_split[i, 'Aa_cells']
          Aa_cells <- stringr::str_split(Aa_cells_character, ",")[[1]]
          aa_cells_character <- metadata_split[i, 'aa_cells']
          aa_cells <- stringr::str_split(aa_cells_character, ",")[[1]]
          genes <- metadata_split[i, 'GeneList']
          genes <- stringr::str_split(genes, ",")[[1]]
          gene.cnt <- 0

          #snp
          results_snp <- data.frame(SNVid = character(),
                                    group = character(),
                                    Geneid = character(),
                                    sample_size_1 = integer(),
                                    sample_size_2 = integer(),
                                    sample_size_3 = integer(),
                                    theta_1 = double(),
                                    theta_2 = double(),
                                    theta_3 = double(),
                                    mu_1 = double(),
                                    mu_2 = double(),
                                    mu_3 = double(),
                                    size_1 = double(),
                                    size_2 = double(),
                                    size_3 = double(),
                                    prob_1 = double(),
                                    prob_2 = double(),
                                    prob_3 = double(),
                                    total_mean_1 = double(),
                                    total_mean_2 = double(),
                                    total_mean_3 = double(),
                                    chi = double(),
                                    pvalue = double(),
                                    adjusted_pvalue = double(),
                                    Remark = character(),
                                    stringsAsFactors=FALSE)

          #gene
          for (gene in genes){
            gene.cnt <- gene.cnt + 1
            counts_1 <- unlist(expressionMatrix[gene, AA_cells]) # gene expression for AA group
            counts_2 <- unlist(expressionMatrix[gene, Aa_cells]) # gene expression for Aa group
            counts_3 <- unlist(expressionMatrix[gene, aa_cells]) # gene expression for aa group
            results_gene <- data.frame(group = j, SNVid = snvid, Geneid = gene, sample_size_1 = length(counts_1), sample_size_2 = length(counts_2), sample_size_3 = length(counts_3),
                                       theta_1 = NA, theta_2 = NA, theta_3 = NA, mu_1 = NA, mu_2 = NA, mu_3 = NA, size_1 = NA, size_2 = NA, size_3 = NA, prob_1 = NA, prob_2 = NA, prob_3 = NA,
                                       total_mean_1 = NA, total_mean_2 = NA, total_mean_3 = NA, chi = NA, pvalue = NA, adjusted_pvalue = NA, Remark = NA,
                                       stringsAsFactors=FALSE)

            #calculate mean
            totalMean_1 <- mean(counts_1)
            totalMean_2 <- mean(counts_2)
            totalMean_3 <- mean(counts_3)

            #filling data
            results_gene[1,"total_mean_1"] <- totalMean_1
            results_gene[1,"total_mean_2"] <- totalMean_2
            results_gene[1,"total_mean_3"] <- totalMean_3

            #build ZINB model
            zinb <- function(counts) {
              # 如果表达数据全都是0值,对ZINB模型的参数直接赋值
              if (sum(counts == 0) == length(counts)) {
                theta <- 1
                mu <- 0
                size <- 1
                prob <- size / (size + mu)
              } else {
                # 如果表达数据不全是0值，第一次建模，用gamlssML
                zinb_gamlssML <- try(gamlssML(counts, family = "ZINBI"), silent = TRUE)
                if ('try-error' %in% class(zinb_gamlssML)) {
                  # 第一次建模失败，第二次建模用zeroinfl
                  zinb_zeroinfl <- try(zeroinfl(formula = counts ~ 1 | 1, dist = "negbin"), silent = TRUE)
                  if ('try-error' %in% class(zinb_zeroinfl)) {
                    # 如果第二次建模还失败，就在results_gene的Remark列标注
                    print("MLE of ZINB failed!")
                    return(list(theta = NA, mu = NA, size = NA, prob = NA))
                  } else {
                    # 如果第二次建模没有失败，就从建模结果中得到参数
                    zinb <- zinb_zeroinfl
                    theta <- plogis(zinb$coefficients$zero)
                    mu <- exp(zinb$coefficients$count)
                    size <- zinb$theta
                    prob <- size / (size + mu)
                  }
                } else {
                  # 如果第一次建模没有失败，就从建模结果中得到参数
                  zinb <- zinb_gamlssML
                  theta <- zinb$nu
                  mu <- zinb$mu
                  size <- 1 / zinb$sigma
                  prob <- size / (size + mu)
                }
              }
              return(list(theta = theta, mu = mu, size = size, prob = prob))
            }

            zinb(counts_1) # 调用qtl-calling函数，传参
            zinb(counts_2)
            zinb(counts_3)

            #estimate AA params
            params_1 <- zinb(counts_1)
            theta_1 <- params_1[['theta']]
            mu_1 <- params_1[['mu']]
            size_1 <- params_1[['size']]
            prob_1 <- params_1[['prob']]
            #print(paste("counts_1",counts_1,"theta_1:", theta_1, "mu_1:", mu_1, "size_1:", size_1, "prob_1:", prob_1))

            #estimate Aa params
            params_2 <- zinb(counts_2)
            theta_2 <- params_2[['theta']]
            mu_2 <- params_2[['mu']]
            size_2 <- params_2[['size']]
            prob_2 <- params_2[['prob']]
            #print(paste("counts_2",counts_2,"theta_2:", theta_2, "mu_2:", mu_2, "size_2:", size_2, "prob_2:", prob_2))

            #estimate aa params
            params_3 <- zinb(counts_3)
            theta_3 <- params_3[['theta']]
            mu_3 <- params_3[['mu']]
            size_3 <- params_3[['size']]
            prob_3 <- params_3[['prob']]
            #print(paste("counts_3",counts_3,"theta_3:", theta_3, "mu_3:", mu_3, "size_3:", size_3, "prob_3:", prob_3))

            #combine AA,Aa and aa
            params_combined <- zinb(c(counts_1, counts_2, counts_3))
            theta_res <- params_combined[['theta']]
            mu_res <- params_combined[['mu']]
            size_res <- params_combined[['size']]
            prob_res <- params_combined[['prob']]
            #print(paste("theta_res:", theta_res, "mu_res:", mu_res, "size_res:", size_res, "prob_res:", prob_res))

            #calculate p-value
            logL <- function(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2, counts_3, theta_3, size_3, prob_3){
              logL_1 <- sum(dzinegbin(counts_1, size = size_1, prob = prob_1, pstr0 = theta_1, log = TRUE)) # log-likelihood for count1 under parameter
              logL_2 <- sum(dzinegbin(counts_2, size = size_2, prob = prob_2, pstr0 = theta_2, log = TRUE)) # log-likelihood for count2 under parameter
              logL_3 <- sum(dzinegbin(counts_3, size = size_3, prob = prob_3, pstr0 = theta_3, log = TRUE)) # log-likelihood for count3 under parameter
              logL <- logL_1 + logL_2 + logL_3
              logL
            }

            logL_1 <- logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2, size_2, prob_2, counts_3, theta_3, size_3, prob_3)
            logL_2 <- logL(counts_1, theta_res, size_res, prob_res, counts_2, theta_res, size_res, prob_res, counts_3, theta_res, size_res, prob_res)
            chi <- logL_1 - logL_2

            pvalue <- try(1 - pchisq(2 * chi , df = 4))
            if (class(pvalue) == "try-error") return(NA)

            #filling data into data frame
            results_gene[1,"theta_1"] <- theta_1
            results_gene[1,"theta_2"] <- theta_2
            results_gene[1,"theta_3"] <- theta_3
            results_gene[1,"mu_1"] <- mu_1
            results_gene[1,"mu_2"] <- mu_2
            results_gene[1,"mu_3"] <- mu_3
            results_gene[1,"size_1"] <- size_1
            results_gene[1,"size_2"] <- size_2
            results_gene[1,"size_3"] <- size_3
            results_gene[1,"prob_1"] <- prob_1
            results_gene[1,"prob_2"] <- prob_2
            results_gene[1,"prob_3"] <- prob_3
            results_gene[1,"chi"] <- chi
            results_gene[1,"pvalue"] <- pvalue

            results_snp <- rbind(results_snp, results_gene)
          }

          #return
          return(results_snp)

          #eQTLcalling end here
        }

        #final result
        result <- data.frame(
          group = character(),
          SNVid = character(),
          Geneid = character(),
          sample_size_1 = integer(),
          sample_size_2 = integer(),
          sample_size_3 = integer(),
          theta_1 = double(),
          theta_2 = double(),
          theta_3 = double(),
          mu_1 = double(),
          mu_2 = double(),
          mu_3 = double(),
          size_1 = double(),
          size_2 = double(),
          size_3 = double(),
          prob_1 = double(),
          prob_2 = double(),
          prob_3 = double(),
          total_mean_1 = double(),
          total_mean_2 = double(),
          total_mean_3 = double(),
          chi = double(),
          pvalue = double(),
          adjusted_pvalue = double(),
          Remark = character(),
          stringsAsFactors=FALSE)

        #process visible
        message("start calling eQTL")

        metadata_split1 <- metadata$group == j
        metadata_split <- metadata[metadata_split1, ]

        for(i in 1:dim(metadata_split)[1]){
          message(paste0("processing ", (i / dim(metadata_split)[1])*100,"%"))
          result <- rbind(result, eQTLcalling(i))
        }
        message("finished!")

        #change column names
        colnames(result)[colnames(result) == 'sample_size_1'] <- 'sample_size_AA'
        colnames(result)[colnames(result) == 'sample_size_2'] <- 'sample_size_Aa'
        colnames(result)[colnames(result) == 'sample_size_3'] <- 'sample_size_aa'
        colnames(result)[colnames(result) == 'theta_1'] <- 'theta_AA'
        colnames(result)[colnames(result) == 'theta_2'] <- 'theta_Aa'
        colnames(result)[colnames(result) == 'theta_3'] <- 'theta_aa'
        colnames(result)[colnames(result) == 'mu_1'] <- 'mu_AA'
        colnames(result)[colnames(result) == 'mu_2'] <- 'mu_Aa'
        colnames(result)[colnames(result) == 'mu_3'] <- 'mu_aa'
        colnames(result)[colnames(result) == 'size_1'] <- 'size_AA'
        colnames(result)[colnames(result) == 'size_2'] <- 'size_Aa'
        colnames(result)[colnames(result) == 'size_3'] <- 'size_aa'
        colnames(result)[colnames(result) == 'prob_1'] <- 'prob_AA'
        colnames(result)[colnames(result) == 'prob_2'] <- 'prob_Aa'
        colnames(result)[colnames(result) == 'prob_3'] <- 'prob_aa'
        colnames(result)[colnames(result) == 'total_mean_1'] <- 'total_mean_AA'
        colnames(result)[colnames(result) == 'total_mean_2'] <- 'total_mean_Aa'
        colnames(result)[colnames(result) == 'total_mean_3'] <- 'total_mean_aa'
      }else{
        stop("sparcity can only be selected as 'TRUE' or 'FALSE'")
      }


      # adjust p-value
      result[,"adjusted_pvalue"] <- p.adjust(result[,"pvalue"], method = p.adjust.method)
      # order adjust p-value
      result <- result[order(result[,"adjusted_pvalue"]),]
      # 重新编号
      rownames(result) <- NULL
      # 根据阈值过滤SNV-gene pairs
      result <- result[result$adjusted_pvalue <= p.adjust.Threshold, ]
      result_all <- rbind(result_all, result)
    }
    return(result_all)

  }


  POSSION <- function(expressionMatrix, snvMatrix = NULL, metadata, sparcity = FALSE, p.adjust.method = "bonferroni", p.adjust.Threshold = 0.05){

    unique_group <- unique(metadata$group)

    result_all <- data.frame()

    for(j in unique_group){
      if(sparcity == TRUE){
        eQTLcalling <- function(i){
          if (i %% 100 == 0) {
            gc()
          }

          # get data from metadata : rsID,snvid,cells,ref_cells,alt_cells,genes ------
          # set gene.cnt = 0
          snvid <- metadata_split[i, "SNVid"]
          cells_character <- metadata_split[i, 'CellList']
          cells <- stringr::str_split(cells_character, ",")[[1]]
          ref_cells_character <- metadata_split[i, 'Ref_cells']
          ref_cells <- stringr::str_split(ref_cells_character, ",")[[1]]
          alt_cells_character <- metadata_split[i, 'Alt_cells']
          alt_cells <- stringr::str_split(alt_cells_character, ",")[[1]]
          genes <- metadata_split[i, 'GeneList']
          genes <- stringr::str_split(genes, ",")[[1]]
          gene.cnt <- 0

          # snp ----------------------------------------------------------------------
          results_snp <- data.frame(
            SNVid = character(),
            Geneid = character(),
            sample_size_1 = integer(),
            sample_size_2 = integer(),
            lambda_1 = double(),
            lambda_2 = double(),
            total_mean_1 = double(),
            total_mean_2 = double(),
            foldChange = double(),
            chi_pos = double(),
            pvalue_pos = double(),
            adjusted_pvalue_pos = double(),
            Remark = character(),
            stringsAsFactors = FALSE)

          # gene ---------------------------------------------------------------------
          for (gene in genes){
            gene.cnt <- gene.cnt + 1
            # gene expression for ref group
            counts_1 <- unlist(expressionMatrix[gene, ref_cells])  # 从表达矩阵中选择指定基因在ref_cells中的counts值
            # gene expression for alt group
            counts_2 <- unlist(expressionMatrix[gene, alt_cells])
            results_gene <- data.frame(
              group = j,
              SNVid = snvid,
              Geneid = gene,
              sample_size_1 = length(counts_1),
              sample_size_2 = length(counts_2),
              lambda_1 = NA,
              lambda_2 = NA,
              total_mean_1 = NA,
              total_mean_2 = NA,
              foldChange = NA,
              chi_pos = NA,
              pvalue_pos = NA,
              adjusted_pvalue_pos = NA,
              Remark = NA,
              stringsAsFactors = FALSE)

            # calculate fold change --------------------------------------------------
            totalMean_1 <- mean(counts_1)
            totalMean_2 <- mean(counts_2)
            foldChange  <- totalMean_1/totalMean_2

            # filling data -----------------------------------------------------------
            results_gene[1,"total_mean_1"] <- totalMean_1
            results_gene[1,"total_mean_2"] <- totalMean_2
            results_gene[1,"foldChange"]   <- foldChange


            # build  model -------------------------------------------------------

            # 这里要在我们的说明文档说明，gamlssML、zeroinfl和glmmTMB这三个函数适用于有零值存在的情况，
            # 而后面三种函数适用于没有零值存在的情况。提醒用户根据数据内容选择。

            build_model <- function(counts) {

              options(show.error.messages = FALSE)  # 在代码运行时不显示错误消息
              possion_glm <- try(glm(formula = counts ~ 1,
                                     family = "poisson",
                                     data= data.frame(counts)), silent = TRUE)
              options(show.error.messages = TRUE)  # 在运行代码时显示错误消息

              if('try-error' %in% class(possion_glm)){
                print("MLE of possion failed! Please choose another model")
                return(lambda = NA)

              }else{
                # 拟合成功，获取参数
                poisson <- possion_glm
                lambda <- poisson$coefficients
              }

              return(list(lambda = lambda))
            }

            # estimate Ref params ----------------------------------------------------
            params_1 <- build_model(counts_1) # 调用qtl-calling函数，传参
            lambda_1 <- params_1[['lambda']]
            # print(paste("counts_1",counts_1,"theta_1:", theta_1, "mu_1:", mu_1, "size_1:", size_1, "prob_1:", prob_1))


            #estimate Alt params -----------------------------------------------------
            params_2 <- build_model(counts_2)
            lambda_2 <- params_2[['lambda']]
            #print(paste("counts_2",counts_2,"theta_2:", theta_2, "mu_2:", mu_2, "size_2:", size_2, "prob_2:", prob_2))

            #combine Ref and Alt，得到综合的模型参数 -----------------------------------------------------
            params_combined <- build_model(c(counts_1, counts_2))
            lambda_res <- params_combined[['lambda']]
            #print(paste("theta_res:", theta_res, "mu_res:", mu_res, "size_res:", size_res, "prob_res:", prob_res))


            # possion
            # 定义函数logL，计算给定参数下的对数似然函数值
            logL_pos <- function(counts_a, lambda1, counts_b, lambda2){
              # 计算两个模型的对数似然函数值
              logL_A <- sum(dpois(counts_a, lambda1, log = TRUE))
              logL_B <- sum(dpois(counts_b, lambda2, log = TRUE))
              logL <- logL_A + logL_B
              logL
            }

            # 用上面定义的logL函数计算ref模型和alt模型的对数似然函数值
            logL_pos1 <- logL_pos(counts_1, lambda_1, counts_2, lambda_2)
            # 用上面定义的logL函数计算counts_1模型和counts_2模型用综合参数的对数似然函数值
            logL_pos2 <- logL_pos(counts_1, lambda_res, counts_2, lambda_res)

            # 得到卡方统计量 chi，这个统计量用于衡量两个模型的拟合优度，进而检验它们是否有显著差异
            chi_pos <- logL_pos1 - logL_pos2
            # 通过 pchisq 函数，基于卡方统计量和自由度，得到p值，这个p值表示了在两个分布模型之间的对数似然比差异的显著性
            pvalue_pos <- 1 - pchisq(2 * chi_pos , df = 3)



            # filling data into data frame ------------------------------

            results_gene[1,"lambda_1"] <- lambda_1
            results_gene[1,"lambda_2"] <- lambda_2

            results_gene[1,"chi_pos"] <- chi_pos
            results_gene[1,"pvalue_pos"] <- pvalue_pos

            results_snp <- rbind(results_snp, results_gene)
          }

          # return
          return(results_snp)

          # eQTLcalling end here
        }

        # final result
        result <- data.frame(
          SNVid = character(),
          Geneid = character(),
          sample_size_1 = integer(),
          sample_size_2 = integer(),
          lambda_1 = double(),
          lambda_2 = double(),
          total_mean_1 = double(),
          total_mean_2 = double(),
          foldChange = double(),
          chi_pos = double(),
          pvalue_pos = double(),
          adjusted_pvalue = double(),
          Remark = character(),
          stringsAsFactors=FALSE)

        #process visible
        message("start calling eQTL")

        metadata_split1 <- metadata$group == j
        metadata_split <- metadata[metadata_split1, ]

        for(i in 1:dim(metadata_split)[1]){
          message(paste0("processing ", (i / dim(metadata_split)[1])*100,"%"))
          # 调用eQTLcalling函数传参给eQTLsc function
          result <- rbind(result, eQTLcalling(i))
        }
        message("finished!")


        #change column names
        colnames(result)[colnames(result) == 'sample_size_1'] <- 'sample_size_Ref'
        colnames(result)[colnames(result) == 'sample_size_2'] <- 'sample_size_Alt'

        colnames(result)[colnames(result) == 'lambda_1'] <- 'lambda_Ref'
        colnames(result)[colnames(result) == 'lambda_2'] <- 'lambda_Alt'
        colnames(result)[colnames(result) == 'total_mean_1'] <- 'total_mean_Ref'
        colnames(result)[colnames(result) == 'total_mean_2'] <- 'total_mean_Alt'
      }else if(sparcity == FALSE){
        eQTLcalling <- function(i){
          if (i %% 100 == 0) {
            gc()
          }

          # get data from metadata : rsID,snvid,cells,ref_cells,alt_cells,genes ------
          # set gene.cnt = 0
          snvid <- metadata_split[i, "SNVid"]
          cells_character <- metadata_split[i, 'CellList']
          cells <- stringr::str_split(cells_character, ",")[[1]]
          AA_cells_character <- metadata_split[i, 'AA_cells']
          AA_cells <- stringr::str_split(AA_cells_character, ",")[[1]]
          Aa_cells_character <- metadata_split[i, 'Aa_cells']
          Aa_cells <- stringr::str_split(Aa_cells_character, ",")[[1]]
          aa_cells_character <- metadata_split[i, 'aa_cells']
          aa_cells <- stringr::str_split(aa_cells_character, ",")[[1]]
          genes <- metadata_split[i, 'GeneList']
          genes <- stringr::str_split(genes, ",")[[1]]
          gene.cnt <- 0

          # snp ----------------------------------------------------------------------
          results_snp <- data.frame(SNVid = character(),
                                    group = character(),
                                    Geneid = character(),
                                    sample_size_1 = integer(),
                                    sample_size_2 = integer(),
                                    sample_size_3 = integer(),
                                    lambda_1 = double(),
                                    lambda_2 = double(),
                                    lambda_3 = double(),
                                    total_mean_1 = double(),
                                    total_mean_2 = double(),
                                    total_mean_3 = double(),
                                    chi_pos = double(),
                                    pvalue_pos = double(),
                                    adjusted_pvalue_pos = double(),
                                    Remark = character(),
                                    stringsAsFactors=FALSE)

          # gene ---------------------------------------------------------------------
          for (gene in genes){
            gene.cnt <- gene.cnt + 1
            counts_1 <- unlist(expressionMatrix[gene, AA_cells]) # gene expression for AA group
            counts_2 <- unlist(expressionMatrix[gene, Aa_cells]) # gene expression for Aa group
            counts_3 <- unlist(expressionMatrix[gene, aa_cells]) # gene expression for aa group
            results_gene <- data.frame(
              group = j,
              SNVid = snvid,
              Geneid = gene,
              sample_size_1 = length(counts_1),
              sample_size_2 = length(counts_2),
              sample_size_3 = length(counts_3),
              lambda_1 = NA,
              lambda_2 = NA,
              lambda_3 = NA,
              total_mean_1 = NA,
              total_mean_2 = NA,
              total_mean_3 = NA,
              chi_pos = NA,
              pvalue_pos = NA,
              adjusted_pvalue_pos = NA,
              Remark = NA,
              stringsAsFactors=FALSE)

            # calculate fold change --------------------------------------------------
            totalMean_1 <- mean(counts_1)
            totalMean_2 <- mean(counts_2)
            totalMean_3 <- mean(counts_3)

            # filling data -----------------------------------------------------------
            results_gene[1,"total_mean_1"] <- totalMean_1
            results_gene[1,"total_mean_2"] <- totalMean_2
            results_gene[1,"total_mean_3"] <- totalMean_3


            # build  model -------------------------------------------------------

            # 这里要在我们的说明文档说明，gamlssML、zeroinfl和glmmTMB这三个函数适用于有零值存在的情况，
            # 而后面三种函数适用于没有零值存在的情况。提醒用户根据数据内容选择。

            build_model <- function(counts) {

              options(show.error.messages = FALSE)  # 在代码运行时不显示错误消息
              possion_glm <- try(glm(formula = counts ~ 1,
                                     family = "poisson",
                                     data= data.frame(counts)), silent = TRUE)
              options(show.error.messages = TRUE)  # 在运行代码时显示错误消息

              if('try-error' %in% class(possion_glm)){
                print("MLE of possion failed! Please choose another model")
                return(lambda = NA)

              }else{
                # 拟合成功，获取参数
                poisson <- possion_glm
                lambda <- poisson$coefficients
              }

              return(list(lambda = lambda))
            }

            # estimate AA params ----------------------------------------------------
            params_1 <- build_model(counts_1) # 调用qtl-calling函数，传参
            lambda_1 <- params_1[['lambda']]
            # print(paste("counts_1",counts_1,"theta_1:", theta_1, "mu_1:", mu_1, "size_1:", size_1, "prob_1:", prob_1))


            #estimate Aa params -----------------------------------------------------
            params_2 <- build_model(counts_2)
            lambda_2 <- params_2[['lambda']]
            #print(paste("counts_2",counts_2,"theta_2:", theta_2, "mu_2:", mu_2, "size_2:", size_2, "prob_2:", prob_2))

            # estimate aa params ----------------------------------------------------
            params_3 <- build_model(counts_3)
            lambda_3 <- params_3[['lambda']]

            #得到综合的模型参数 -----------------------------------------------------
            params_combined <- build_model(c(counts_1, counts_2, counts_3))
            lambda_res <- params_combined[['lambda']]
            #print(paste("theta_res:", theta_res, "mu_res:", mu_res, "size_res:", size_res, "prob_res:", prob_res))


            # possion
            # 定义函数logL，计算给定参数下的对数似然函数值
            logL_pos <- function(counts_a, lambda1, counts_b, lambda2, counts_c, lambda3){
              # 计算两个模型的对数似然函数值
              logL_A <- sum(dpois(counts_a, lambda1, log = TRUE))
              logL_B <- sum(dpois(counts_b, lambda2, log = TRUE))
              logL_C <- sum(dpois(counts_c, lambda3, log = TRUE))
              logL <- logL_A + logL_B + logL_C
              logL
            }

            # 用上面定义的logL函数计算ref模型和alt模型的对数似然函数值
            logL_pos1 <- logL_pos(counts_1, lambda_1, counts_2, lambda_2, counts_3, lambda_3)
            # 用上面定义的logL函数计算counts_1模型和counts_2模型用综合参数的对数似然函数值
            logL_pos2 <- logL_pos(counts_1, lambda_res, counts_2, lambda_res, counts_3, lambda_res)

            # 得到卡方统计量 chi，这个统计量用于衡量两个模型的拟合优度，进而检验它们是否有显著差异
            chi_pos <- logL_pos1 - logL_pos2
            # 通过 pchisq 函数，基于卡方统计量和自由度，得到p值，这个p值表示了在两个分布模型之间的对数似然比差异的显著性
            pvalue_pos <- 1 - pchisq(2 * chi_pos , df = 3)



            # filling data into data frame ------------------------------

            results_gene[1,"lambda_1"] <- lambda_1
            results_gene[1,"lambda_2"] <- lambda_2
            results_gene[1,"lambda_3"] <- lambda_3

            results_gene[1,"chi_pos"] <- chi_pos
            results_gene[1,"pvalue_pos"] <- pvalue_pos

            results_snp <- rbind(results_snp, results_gene)
          }

          # return
          return(results_snp)

          # eQTLcalling end here
        }

        # final result
        result <- data.frame(SNVid = character(),
                             group = character(),
                             Geneid = character(),
                             sample_size_1 = integer(),
                             sample_size_2 = integer(),
                             sample_size_3 = integer(),
                             lambda_1 = double(),
                             lambda_2 = double(),
                             lambda_3 = double(),
                             total_mean_1 = double(),
                             total_mean_2 = double(),
                             total_mean_3 = double(),
                             chi_pos = double(),
                             pvalue_pos = double(),
                             adjusted_pvalue_pos = double(),
                             Remark = character(),
                             stringsAsFactors=FALSE)
        #process visible
        message("start calling eQTL")

        metadata_split1 <- metadata$group == j
        metadata_split <- metadata[metadata_split1, ]

        for(i in 1:dim(metadata_split)[1]){
          message(paste0("processing ", (i / dim(metadata_split)[1])*100,"%"))
          # 调用eQTLcalling函数传参给eQTLsc function
          result <- rbind(result, eQTLcalling(i))
        }
        message("finished!")


        #change column names
        colnames(result)[colnames(result) == 'sample_size_1'] <- 'sample_size_AA'
        colnames(result)[colnames(result) == 'sample_size_2'] <- 'sample_size_Aa'
        colnames(result)[colnames(result) == 'sample_size_3'] <- 'sample_size_aa'

        colnames(result)[colnames(result) == 'lambda_1'] <- 'lambda_AA'
        colnames(result)[colnames(result) == 'lambda_2'] <- 'lambda_Aa'
        colnames(result)[colnames(result) == 'lambda_3'] <- 'lambda_aa'
        colnames(result)[colnames(result) == 'total_mean_1'] <- 'total_mean_AA'
        colnames(result)[colnames(result) == 'total_mean_2'] <- 'total_mean_Aa'
        colnames(result)[colnames(result) == 'total_mean_3'] <- 'total_mean_aa'
      }else{
        stop("请输入正确的snvtype")
      }

      # 加p值矫正方法
      if (!p.adjust.method %in% c("bonferroni", "holm", "hochberg", "hommel", "BH")) {
        stop("Invalid p-adjusted method.
           Please choose from 'bonferroni', 'holm', 'hochberg', 'hommel', or'fdr or BH'.")
      }

      # adjust p-value
      result[,"adjusted_pvalue_pos"] <- p.adjust(result[,"pvalue_pos"], method = p.adjust.method)
      # order adjust p-value
      result <- result[order(result[,"adjusted_pvalue_pos"]),]
      # 重新编号
      rownames(result) <- NULL
      # 根据阈值过滤SNV-gene pairs
      result <- result[result$adjusted_pvalue_pos <= p.adjust.Threshold, ]

      if(exists("lastFuncGrad") & exists("lastFuncParam")){
        lastFuncGrad <- NULL
        lastFuncParam <- NULL
        remove(lastFuncGrad, lastFuncParam, envir=.GlobalEnv)
      }
      result_all <- rbind(result_all, result)
    }
    return(result_all)
  }


  if(useModel == "zinb"){
    result <- ZINB(expressionMatrix, metadata = metadata, sparcity = sparcity, p.adjust.method = p.adjust.method, p.adjust.Threshold = p.adjust.Threshold)
    return(result)
  }else if(useModel == "possion"){
    result <- POSSION(expressionMatrix, metadata = metadata, sparcity = sparcity, p.adjust.method = p.adjust.method, p.adjust.Threshold = p.adjust.Threshold)
    return(result)
  }else if(useModel == "linear"){
    result <- LINEAR(expressionMatrix, snvMatrix, metadata, sparcity, p.adjust.method = p.adjust.method, p.adjust.Threshold = p.adjust.Threshold)
    return(result)
  }else{
    stop("Invalid model Please choose from 'zinb', 'possion' , or 'linear'.")
  }
}
