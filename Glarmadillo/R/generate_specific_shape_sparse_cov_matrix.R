#' @title Generate Specific Shape Sparse Covariance Matrix
#' @description
#' Generates a covariance matrix and corresponding data matrix (Y)
#' with a specific shape defined by a given shape matrix (M).
#' This function is particularly useful for simulating data with
#' predefined covariance structures, facilitating the testing of
#' statistical methods such as sparse covariance estimation.
#'
#' @param n The number of variables (rows of Y and dimensions of M).
#' @param p The number of samples (columns of Y).
#' @param M The shape matrix used to define the structure of the covariance matrix.
#'          Must be a positive definite square matrix of size n x n.
#'
#' @return A list containing two elements:
#'         - Y: A n by p data matrix, where each column represents a sample,
#'              and each row represents a variable.
#'         - cov_Y: The n by n covariance matrix of the transposed data matrix Y.
#'          This covariance matrix reflects the structure imposed by the shape matrix M.
#' @importFrom stats rnorm
#' @examples
#' # Generate a 10x10 specific shape sparse covariance matrix
#' shape_matrix <- matrix(rnorm(100), 10, 10)
#' shape_matrix <- shape_matrix %*% t(shape_matrix) # Making it positive definite
#' result <- generate_specific_shape_sparse_cov_matrix(n = 10, p = 5, M = shape_matrix)
#' Y <- result$Y
#' cov_Y <- result$cov_Y
#'
#' @export
generate_specific_shape_sparse_cov_matrix <- function(n, p, M) {

  M <- as.matrix(M)

  if (!is.matrix(M) || nrow(M) != ncol(M) || nrow(M) != n) {
    stop("M must be a square matrix of size n x n")
  }

  if (p > n) {
    stop("p must not be greater than n")
  }

  # 检查 M 是否可进行 Cholesky 分解
  if (!isSymmetric(M) || any(eigen(M)$values <= 0)) {
    stop("M must be a positive definite matrix for Cholesky decomposition")
  }

  # 使用M进行Cholesky分解
  L <- tryCatch({
    chol(M)
  }, error = function(e) {
    stop("Cholesky decomposition failed: ", e$message)
  })

  # 生成随机数据矩阵 Z
  Z <- replicate(p, rnorm(n))

  # 计算 Y = solve(L, Z)
  Y <- matrix(solve(L, Z),n,p)

  # 计算协方差矩阵 cov(Y)
  cov_Y <- cov(t(Y))

  return(list(Y = Y, cov_Y = cov_Y))
}
