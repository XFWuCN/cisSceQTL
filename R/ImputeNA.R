#' Title
#'
#' @param expressionMatrix 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
ImputeNA <- function(expressionMatrix, method = "linear"){
  if (sum(is.na(expressionMatrix)) > 0) {
    if(method == "linear"){
      # 使用线性插值填充缺失值
      expressionMatrix.linear <- expressionMatrix
      expressionMatrix.linear[is.na(expressionMatrix.linear)] <- approx(1:ncol(expressionMatrix.linear), 1:nrow(expressionMatrix.linear),
                                                                        xout = which(is.na(expressionMatrix.linear)),
                                                                        y = expressionMatrix.linear[!is.na(expressionMatrix.linear)],
                                                                        method = "linear", rule = 2)$y
      return(expressionMatrix.linear)
    } else if(method == "polynomial"){
      # 使用多项式插值填充缺失值
      expressionMatrix.polynomial <- expressionMatrix
      expressionMatrix.polynomial[is.na(expressionMatrix.polynomial)] <- spline(1:ncol(expressionMatrix.polynomial),
                                                                                expressionMatrix.polynomial[!is.na(expressionMatrix.polynomial)],
                                                                                n = nrow(expressionMatrix.polynomial))$y
      return(expressionMatrix.polynomial)
    } else if(method == "knn"){
      # 使用K近邻插值填充缺失值
      expressionMatrix.knn <- expressionMatrix
      library(impute)
      expressionMatrix.knn <- na.knn(expressionMatrix.knn, k = 5)
      return(expressionMatrix.knn)
    } else if(method == "spline"){
      # 使用样条插值填充缺失值
      expressionMatrix.spline <- expressionMatrix
      expressionMatrix.spline[is.na(expressionMatrix.spline)] <- spline(1:ncol(expressionMatrix.spline),
                                                                        expressionMatrix.spline[!is.na(expressionMatrix.spline)],
                                                                        n = nrow(expressionMatrix.spline))$y
      return(expressionMatrix.spline)
    } else {
      stop("Invalid method. Please choose from 'linear', 'polynomial', 'knn', or 'spline'.")
    }
  } else {
    cat("Your data does not contain any missing values!\n")
    return(expressionMatrix)
  }
}
