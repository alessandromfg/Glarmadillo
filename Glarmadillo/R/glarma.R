# glarma.R
# Part of the Glarmadillo package
# Copyright (C) 2023 [Alessandro]
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#' @title Solve Graphical Lasso with Armadillo
#' @description
#' This function solves the Graphical Lasso (GLasso) problem using the Armadillo library.
#' GLasso is a technique used in statistical learning and network analysis to
#' estimate sparse inverse covariance matrices from observed data.
#'
#' @param s A symmetric, positive-definite sample covariance matrix.
#'          It should be a square matrix representing the covariance matrix of the variables.
#' @param rho A positive scalar representing the regularization parameter.
#'            It controls the sparsity level of the inverse covariance matrix.
#' @param mtol A numeric value representing the convergence threshold for the main algorithm.
#'             It determines the condition under which the iterative process will stop.
#'             Default is 1e-4.
#' @param maxIterations An integer value specifying the maximum number of iterations
#'                      allowed for the algorithm. Default is 10000.
#' @param ltol A numeric value representing the convergence threshold for the Lasso solver.
#'             It is used to control the Lasso solving process within the algorithm.
#'             Default is 1e-6.
#'
#' @return Returns a sparse inverse covariance matrix estimated by solving
#'         the Graphical Lasso problem. The sparsity is controlled by the 'rho' parameter.
#'
#' @export
#' @importFrom stats runif
#' @examples
#' # Generate a sample covariance matrix
#' s <- matrix(runif(100), nrow = 10)
#' s <- t(s) %*% s
#' # Solve the Graphical Lasso problem with default parameters
#' inv_cov_matrix <- glarma(s, rho = 0.1)
#' # Solve with custom convergence thresholds and maximum iterations
#' inv_cov_matrix <- glarma(s, rho = 0.1, mtol = 1e-5, maxIterations = 5000, ltol = 1e-6)
#'
#' @useDynLib Glarmadillo
#' @import Rcpp
#' @import RcppArmadillo
glarma <- function(s, rho, mtol = 1e-4, maxIterations = 10000, ltol = 1e-6) {
  # 使用 checkMatrix 函数
  checkMatrix(s)

  # 确保 rho 是非负标量
  if (!is.numeric(rho) || length(rho) != 1 || rho < 0) {
    stop("Regularization parameter 'rho' must be a non-negative scalar.")
  }

  # 如果 rho 为 0，发出警告
  if (rho == 0) {
    warning("With rho=0, there may be convergence problems if the input matrix is not of full rank")
  }

  # 调用 C++ 函数
  result <- glarma_cpp(s, rho, mtol, maxIterations, ltol)

  return(result)
}
#' @title Generate Sparse Covariance Matrix
#' @description
#' Generates a sparse covariance matrix with specified dimension and rank.
#' The generated matrix is scaled to control the magnitude of its elements
#' and can be further sparsified based on a given threshold.
#'
#' @param n The dimension of the covariance matrix (number of rows and columns).
#' @param p The rank of the covariance matrix (number of non-zero eigenvalues).
#'          Must be less than or equal to \code{n}.
#' @param sparse Logical indicating whether to enforce additional sparsity by
#'        setting elements with absolute values below \code{sparse_rho} to zero.
#' @param sparse_rho Numeric threshold for enforcing sparsity.
#'        Only used if \code{sparse} is \code{TRUE}.
#' @param scale_power The exponent used to scale the matrix elements,
#'        with a default value of 1. This helps to control the magnitude
#'        of the covariance matrix elements as \code{n} increases.
#'
#' @return A sparse \code{n} by \code{n} covariance matrix with rank \code{p}.
#'         If \code{sparse} is \code{TRUE}, elements with absolute values below
#'         \code{sparse_rho} are set to zero to increase sparsity,
#'         while ensuring that the matrix is at least semi-definite.
#'
#' @examples
#' # Generate a 10x10 sparse covariance matrix with rank 5
#' sparse_cov_matrix <- generate_sparse_cov_matrix(n = 10, p = 5)
#'
#' # Generate a sparser matrix with elements below 0.3 set to zero
#' sparser_cov_matrix <- generate_sparse_cov_matrix(n = 100, p = 50,
#'                                                  sparse = TRUE, sparse_rho = 0.3)
#'
#' # Generate a matrix with a different scale power
#' sparse_cov_matrix <- generate_sparse_cov_matrix(n = 100, p = 50, scale_power = 2)
#'
#' @export
generate_sparse_cov_matrix <- function(n, p, sparse = TRUE, sparse_rho = 0.3, scale_power = 1) {
  if (p > n) {
    stop("p must not be greater than n")
  }

  # 步骤 1: 在对角线上生成 p 个正随机数和 n-p 个零
  diag_values <- c(runif(p, min = 0.1, max = 1), rep(0, n - p))
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

  # 步骤 5：缩放矩阵元素，使用提供的 scale_power
  S <- S / n^scale_power

  # 如果启用了稀疏性，则进一步将小于sparse_rho的元素置为零
  if (sparse) {
    S[abs(S) < sparse_rho] <- 0
    # 确保矩阵至少是半正定的
    eigen_values <- eigen(S)$values
    min_eigen_value <- min(eigen_values)
    if (min_eigen_value < 0) {
      # 如果有负特征值，进行修正
      S <- S - diag(min_eigen_value, n)
    }
  }

  return(S)
}
#' @title Check Matrix Validity for GLarma
#' @description
#' Internal function to check if a matrix is valid for processing in the glarma function.
#' Specifically, it checks if the matrix is a square, symmetric matrix.
#'
#' @param m A matrix to be checked.
#' @return Returns TRUE if the matrix is valid, otherwise throws an error.
#'
#' @noRd
checkMatrix <- function(m, tolerance = 1e-6) {
  if (!is.matrix(m)) {
    stop("Input is not a matrix.")
  }
  if (nrow(m) != ncol(m)) {
    stop("Matrix is not square.")
  }
  if (!isTRUE(all.equal(m, t(m), tolerance = tolerance))) {
    stop("Matrix is not symmetric.")
  }
  return(TRUE)
}
