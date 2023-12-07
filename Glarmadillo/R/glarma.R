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
#' @return Returns a covariance matrix W and a estimated sparse inverse covariance matrix Theta
#'         estimated by solving the Graphical Lasso problem. The sparsity is controlled
#'         by the 'rho' parameter.
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
