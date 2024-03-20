#' Title
#'
#' @param expressionMatrix A dataframe, describes gene expressions, the row names represent gene IDs and the column names represent cell IDs.
#' @param method Method for imputating NA for gene expression matrix.
#'
#' @return A expressionMatrix with no missing values.
#' @export
#'
#' @examples
#' # Load test data for cisSceQTL
#' data(testGene)
#' # impute the missing values for testGene
#' expressionMatrix <- ImputeNA(testGene, method = "linear")

ImputeNA <- function(expressionMatrix, method = "linear"){
  if (sum(is.na(expressionMatrix)) > 0) {
    if(method == "linear"){
      # impute NA by linear
      expressionMatrix <- expressionMatrix
      expressionMatrix[is.na(expressionMatrix)] <- approx(1:ncol(expressionMatrix.linear), 1:nrow(expressionMatrix.linear),
                                                                        xout = which(is.na(expressionMatrix.linear)),
                                                                        y = expressionMatrix.linear[!is.na(expressionMatrix.linear)],
                                                                        method = "linear", rule = 2)$y
      return(expressionMatrix.linear)
    } else if(method == "polynomial"){
      # impute NA by polynomial
      expressionMatrix.polynomial <- expressionMatrix
      expressionMatrix.polynomial[is.na(expressionMatrix.polynomial)] <- spline(1:ncol(expressionMatrix.polynomial),
                                                                      expressionMatrix.polynomial[!is.na(expressionMatrix.polynomial)],
                                                                                n = nrow(expressionMatrix.polynomial))$y
      return(expressionMatrix.polynomial)
    } else if(method == "knn"){
      # impute NA by knn
      expressionMatrix.knn <- expressionMatrix
      library(impute)
      expressionMatrix.knn <- na.knn(expressionMatrix.knn, k = 5)
      return(expressionMatrix.knn)
    } else if(method == "spline"){
      # impute NA by spline
      expressionMatrix.spline <- expressionMatrix
      expressionMatrix.spline[is.na(expressionMatrix.spline)] <- spline(1:ncol(expressionMatrix.spline),
                                                                        expressionMatrix.spline[!is.na(expressionMatrix.spline)],
                                                                        n = nrow(expressionMatrix.spline))$y
      return(expressionMatrix.spline)
    } else {
      stop("Invalid method. Please choose from 'linear', 'polynomial', 'knn', or 'spline'.")
    }
  } else {
    cat("Your gene expression dataframe does not contain any missing values!\n")
    return(expressionMatrix)
  }
}
