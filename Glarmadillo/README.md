# Glarmadillo

`Glarmadillo` is an R package designed for solving the Graphical Lasso (GLasso) problem using the Armadillo library. It provides an efficient implementation for estimating sparse inverse covariance matrices from observed data, ideal for applications in statistical learning and network analysis.In addition to solving the GLasso problem, this package includes functionality to generate random sparse covariance matrices which are useful in simulations and teaching of graphical models.

## Installation

To install the latest version of `Glarmadillo` from GitHub, run the following commands in R:

```r
# Install the devtools package if it's not already installed
if (!require(devtools)) install.packages("devtools")

# Install the Glarmadillo package from GitHub
devtools::install_github("alessandromfg/Glarmadillo")
```

## Usage
Here's a basic example to demonstrate how to generate a sparse covariance matrix and solve the GLasso problem:

```r
# Load the Glarmadillo package
library(Glarmadillo)

# Set parameters
n <- 20
p <- 5

# Generate a sparse covariance matrix
s1 <- generate_sparse_cov_matrix(n, p, 1, 0.3, 0.5)

# Set the regularization parameter
rho <- 0.1

# Solve the Graphical Lasso problem
result <- glarma(s1, rho, 1e-4, 1000, 1e-7)
```

## Parameter Selection Tips

When selecting parameters for `generate_sparse_cov_matrix`, consider the following guidelines based on the matrix dimension (`n`):

- For smaller dimensions (e.g., `n` around 50), the `scale_power` parameter can be set lower, such as around 0.8.
- For larger dimensions (e.g., `n` equals 1000), `scale_power` can be increased to 1.2. Since the fitting changes with the expansion of `n`, 1.5 is the upper limit and should not be exceeded.
- Regarding sparsity, the sparsity coefficient `sparse_rho` should vary inversely with `scale_power`. For example, if `scale_power` is 0.8, `sparse_rho` can be set to 0.4; if `scale_power` is 1.2, `sparse_rho` can be reduced to 0.1. Adjust according to the specific matrix form.

For `glarma`, the `mtol` parameter, which controls the overall matrix convergence difference, typically does not need adjustment. The number of iterations usually does not reach the maximum, and convergence generally occurs within 20 iterations. The critical aspect is the adjustment of `ltol`. It is recommended to decrease `ltol` as the matrix size increases. For instance, when `n` is 20, `ltol` can be set to 1e-5; when `n` is 1000, it should be set to 1e-7. This is because `ltol` is the convergence condition for each column; if `n` is too large and `ltol` is not sufficiently reduced, the final results can vary significantly.

## License
License This package is free and open-source software licensed under the GNU General Public License version 3 (GPL-3)
