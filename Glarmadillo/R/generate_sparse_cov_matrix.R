#' @title Generate Sparse Covariance Matrix
#' @description
#' Generates a sparse covariance matrix with specified dimension and rank.
#' The generated matrix can be scaled or standardized, and further sparsified
#' based on a given threshold.
#'
#' @param n The dimension of the covariance matrix (number of rows and columns).
#' @param p The rank of the covariance matrix (number of non-zero eigenvalues).
#'          Must be less than or equal to \code{n}.
#' @param standardize Logical indicating whether to standardize the matrix,
#'        setting this to TRUE overrides scale_power and sparse_rho.
#' @param sparse_rho Numeric threshold for enforcing sparsity.
#'        Elements with absolute values below \code{sparse_rho} are set to zero.
#' @param scale_power The exponent used to scale the matrix elements.
#'        Only used if \code{standardize} is \code{FALSE}.
#'
#' @return A \code{n} by \code{n} covariance matrix with rank \code{p}.
#'         If \code{sparse_rho} is greater than zero and \code{standardize} is FALSE,
#'         elements with absolute values below \code{sparse_rho} are set to zero to increase sparsity,
#'         while ensuring that the matrix is at least semi-definite.
#' @importFrom stats cov
#' @examples
#' # Generate a 10x10 sparse covariance matrix with rank 5
#' sparse_cov_matrix <- generate_sparse_cov_matrix(n = 10, p = 5)
#'
#' # Generate a sparser matrix with elements below 0.3 set to zero
#' sparser_cov_matrix <- generate_sparse_cov_matrix(n = 100, p = 50,
#'                                                 sparse_rho = 0.3,
#'                                               standardize = FALSE)
#'
#' # Generate a standardized matrix
#' standardized_cov_matrix <- generate_sparse_cov_matrix(n = 100, p = 50, standardize = TRUE)
#'
#' @export
generate_sparse_cov_matrix <- function(n, p, standardize = TRUE, sparse_rho = 0, scale_power = 0) {
  if (p > n) {
    stop("p must not be greater than n")
  }

  if (standardize && (sparse_rho > 0 || scale_power != 0)) {
    warning("When standardize is TRUE, sparse_rho and scale_power are ignored.")
  }

  # 步骤 1: 在对角线上生成 p 个正随机数和 n-p 个非零的最小值
  min_value <- 1e-6  # 一个小的正数，以避免除以零的情况
  diag_values <- c(runif(p, min = 0.1, max = 1), rep(min_value, n - p))
  diag_values <- sample(diag_values) # 打乱顺序
  A <- diag(diag_values)

  # 步骤 2: 生成一个 n 维随机方阵 B，确保其行列式非零
  B <- matrix(runif(n * n, min = -1, max = 1), n, n)
  while (det(B) == 0) {
    B <- matrix(runif(n * n, min = -1, max = 1), n, n)
  }

  # 步骤 3: 计算正定阵 U = B %*% t(B)
  U <- B %*% t(B)

  # 步骤 4: 计算稀疏协方差矩阵 S = t(U) %*% A %*% U
  S <- t(U) %*% A %*% U

  # 步骤 5：根据选择进行标准化或缩放
  if (standardize) {
    # 生成随机数据矩阵 X
    sqrt_diag_values <- sqrt(diag_values)  # 对角线元素开根号
    X <- U %*% diag(sqrt_diag_values)  # 构造随机数据矩阵 X

    # 标准化数据矩阵 X
    X <- scale(X)  # 默认标准化: 中心化并转换为单位方差

    # 计算标准化后的协方差矩阵，即相关系数矩阵
    S <- cov(X)
  } else {
    if (scale_power != 0) {
      S <- S / n^scale_power
    }
    # 应用稀疏化
    if (sparse_rho > 0) {
      S[abs(S) < sparse_rho] <- 0
      # 保证矩阵是半正定的，将所有负特征值设置为零
      eigen_decomp <- eigen(S)
      eigen_values <- eigen_decomp$values
      eigen_values[eigen_values < 0] <- 0  # 将所有负特征值设置为零
      S <- eigen_decomp$vectors %*% diag(eigen_values) %*% t(eigen_decomp$vectors)
    }
  }

  return(S)
}
