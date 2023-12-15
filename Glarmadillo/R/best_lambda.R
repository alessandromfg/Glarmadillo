#' @title Find Optimal Lambda by Sparsity Level
#' @description
#' This function performs a grid search over a range of lambda values to identify
#' the lambda that achieves a desired level of sparsity in the precision matrix
#' estimated by Graphical Lasso. Sparsity is defined as the proportion of zero
#' elements (excluding the diagonal) in the precision matrix.
#'
#' @param s The sample covariance matrix of the data.
#' @param lambda_grid A numeric vector of lambda values to be tested in the grid search.
#' @param desired_sparsity The target sparsity level as a proportion of zero elements
#'        in the precision matrix. This should be a value between 0 and 1.
#' @param mtol The convergence threshold for Graphical Lasso optimization.
#' @param maxIterations The maximum number of iterations for Graphical Lasso optimization.
#' @param ltol The tolerance for determining whether elements are considered zero
#'        when calculating sparsity.
#'
#' @return A list containing the following components:
#'         - \code{best_lambda}: the lambda value that results in sparsity closest to the desired level.
#'         - \code{best_sparsity_difference}: the smallest difference between achieved and desired sparsity.
#'         - \code{actual_sparsity}: a numeric vector of actual sparsity levels for each lambda tested.
#'         - \code{lambda_grid}: the vector of lambda values tested.
#' @importFrom stats cov
#' @examples
#' # Generate a sparse covariance matrix
#' values <- c(160, 50)
#' n <- values[1]
#' p <- values[2]
#' s <- generate_sparse_cov_matrix(n, p, standardize = TRUE, sparse_rho = 0, scale_power = 0)
#'
#' # Define a sequence of lambda values for the grid search
#' lambda_find <- c(0.1, 0.2, 0.3, 0.4)
#'
#' # Perform a grid search to find the lambda value
#' # that results in a precision matrix with approximately 80% sparsity
#' lambda_results <- find_lambda_by_sparsity(s, lambda_find, desired_sparsity = 0.8)
#'
#' # Inspect the optimal lambda value
#' optimal_lambda <- lambda_results$best_lambda
#'
#' # Inspect the sparsity levels for each lambda tested
#' sparsity_levels <- lambda_results$actual_sparsity
#'
#' @export
find_lambda_by_sparsity <- function(s, lambda_grid, desired_sparsity, mtol = 1e-4, maxIterations = 10000, ltol = 1e-6) {
  best_lambda <- NA
  best_sparsity_difference <- Inf
  sparsity_list <- numeric(length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    result <- glarma(s, lambda, mtol, maxIterations, ltol)
    Theta <- result$Theta

    # 计算非对角线元素中0的比例
    off_diagonal_elements <- Theta[upper.tri(Theta)]
    sparsity <- mean(off_diagonal_elements == 0)
    sparsity_list[i] <- sparsity

    # 计算与期望稀疏度的差异
    sparsity_difference <- abs(sparsity - desired_sparsity)

    if (sparsity_difference < best_sparsity_difference) {
      best_sparsity_difference <- sparsity_difference
      best_lambda <- lambda
    }
  }

  # 返回与期望稀疏度最接近的 lambda 和实际稀疏度
  return(list(
    best_lambda = best_lambda,
    best_sparsity_difference = best_sparsity_difference,
    actual_sparsity = sparsity_list,
    lambda_grid = lambda_grid
  ))
}
