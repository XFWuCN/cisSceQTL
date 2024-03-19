#' Title
#'
#' @param exprsMatrix A dataframe, describes gene expressions, the row names represent gene IDs and the column names represent cell IDs.
#' @param method Method for imputating NA for gene expression dataframe.
#'
#' @return A exprsMatrix with no missing values.
#' @export
#'
#' @examples
#' # Load test data for cisSceQTL
#' data(testGene)
#' # impute the missing values for testGene
#' result <- ImputeNA(ImputeNA, method = "linear")

ImputeNA <- function(exprsMatrix, method = "linear"){
  if (sum(is.na(exprsMatrix)) > 0) {
    if(method == "linear"){
      # impute NA by linear
      exprsMatrix <- exprsMatrix
      exprsMatrix[is.na(exprsMatrix)] <- approx(1:ncol(exprsMatrix.linear), 1:nrow(exprsMatrix.linear),
                                                                        xout = which(is.na(exprsMatrix.linear)),
                                                                        y = exprsMatrix.linear[!is.na(exprsMatrix.linear)],
                                                                        method = "linear", rule = 2)$y
      return(exprsMatrix.linear)
    } else if(method == "polynomial"){
      # impute NA by polynomial
      exprsMatrix.polynomial <- exprsMatrix
      exprsMatrix.polynomial[is.na(exprsMatrix.polynomial)] <- spline(1:ncol(exprsMatrix.polynomial),
                                                                      exprsMatrix.polynomial[!is.na(exprsMatrix.polynomial)],
                                                                                n = nrow(exprsMatrix.polynomial))$y
      return(exprsMatrix.polynomial)
    } else if(method == "knn"){
      # impute NA by knn
      exprsMatrix.knn <- exprsMatrix
      library(impute)
      exprsMatrix.knn <- na.knn(exprsMatrix.knn, k = 5)
      return(exprsMatrix.knn)
    } else if(method == "spline"){
      # impute NA by spline
      exprsMatrix.spline <- exprsMatrix
      exprsMatrix.spline[is.na(exprsMatrix.spline)] <- spline(1:ncol(exprsMatrix.spline),
                                                              exprsMatrix.spline[!is.na(exprsMatrix.spline)],
                                                                        n = nrow(exprsMatrix.spline))$y
      return(exprsMatrix.spline)
    } else {
      stop("Invalid method. Please choose from 'linear', 'polynomial', 'knn', or 'spline'.")
    }
  } else {
    cat("Your gene expression dataframe does not contain any missing values!\n")
    return(exprsMatrix)
  }
}
