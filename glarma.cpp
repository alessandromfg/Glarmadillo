#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Soft-thresholding operator
double softThreshold(double a, double lambda) {
  return std::copysign(std::max(0.0, std::abs(a) - lambda), a);
}

// 初始化W矩阵
arma::mat initializeW(const arma::mat& s, double rho) {
  int n = s.n_rows;
  arma::mat I = arma::eye(n, n);
  arma::mat W = s + rho * I;
  return W;
}

// Lasso求解函数
arma::vec solveLasso(const arma::mat& V, const arma::vec& u, double rho, const arma::mat& s) {
  int p = V.n_rows;
  arma::vec beta = arma::zeros(p);
  arma::vec beta_old(p);
  arma::vec off_diag;

  if (p > 1) {
    off_diag = s.elem(arma::find(arma::trimatu(arma::ones<arma::mat>(p,p), 1)));  // 获取非主对角线的元素
  } else {
    off_diag = arma::zeros<arma::vec>(1); // 如果矩阵大小为1x1，直接设置off_diag为0
  }

  double threshold = 0.001 * arma::mean(arma::abs(off_diag));  // 计算阈值
  double mean_diff;

  do {
    beta_old = beta;

    for (int j = 0; j < p; j++) {
      double numerator = u(j) - arma::dot(V.col(j), beta) + V(j, j) * beta(j);
      beta(j) = softThreshold(numerator, rho) / V(j, j);
    }

    arma::vec diff = beta - beta_old;
    mean_diff = arma::mean(arma::abs(diff));

  } while (mean_diff >= threshold);

  return beta;
}

//处理一元
// Lasso求解函数针对j=1的情况
double solveLassoScalar(double V, double u, double rho, double s12_value) {
  double beta = 0.0;
  double beta_old;
  double threshold = 0.001 * std::abs(s12_value);  // 根据s矩阵的第一行第二列的值来计算阈值
  double diff;

  do {
    beta_old = beta;

    double numerator = u - V * beta;
    beta = softThreshold(numerator, rho) / V;

    diff = std::abs(beta - beta_old);

  } while (diff >= threshold);

  return beta;
}

double maxDifference(const arma::mat& A, const arma::mat& B) {
  return arma::max(arma::max(arma::abs(A - B)));
}


// 修改后的计算\Theta矩阵的函数
arma::mat computeTheta(const arma::mat& W, const arma::mat& Betas) {
  Rcpp::Rcout << "Starting computeTheta...\n";

  int p = W.n_rows;
  arma::mat Theta(p, p, arma::fill::zeros);  // 初始化\Theta为p x p的零矩阵

  for (int j = p-1; j >= 1; j--) {
    arma::vec w12 = W.submat(0, j, j-1, j);
    double w22 = W(j, j);
    arma::vec beta = Betas.col(j).head(j);  // 只取前j个元素作为beta

    double theta22 = 1.0 / (w22 - arma::dot(w12.t(), beta));
    arma::vec theta12 = -beta * theta22;

    Theta.submat(0, j, j-1, j) = theta12;
    Theta(j, j) = theta22;
  }

  Rcpp::Rcout << "Completed computeTheta with Theta:\n" << Theta << "\n";

  return Theta;  // 返回计算得到的\Theta矩阵
}


// 主函数
// [[Rcpp::export]]
List glarma_cpp(const arma::mat& s, double rho) {
  Rcpp::Rcout << "Starting glarma_cpp...\n";

  arma::mat W = initializeW(s, rho);
  arma::mat W_old;

  int p = s.n_rows;
  int maxIterations = 100; // 最大迭代次数
  double tol = 1e-4; // 收敛阈值

  arma::mat Betas(p, p, arma::fill::zeros);  // 初始化Betas矩阵为p x p的零矩阵

  for (int iter = 0; iter < maxIterations; iter++) {
    Rcpp::Rcout << "Starting iteration " << iter << " with W:\n" << W << "\n";

      W_old = W; // 保存W矩阵

      // 处理j=1的情况 - 新的部分
      double W11 = W(0,0);
      double s12_value = s(0,1);
      double beta = solveLassoScalar(W11, s12_value, rho, s12_value);
      Betas(0,1) = beta;  // 保存beta值到Betas矩阵的对应位置
      W(0,1) = W(1,0) = W11 * beta;  // 更新W矩阵的相关元素

      for (int j = p-1; j > 1; j--) {  // 从第p列开始，到第二列结束
        arma::mat W11 = W.submat(0, 0, j-1, j-1); // 获取W子矩阵
        arma::vec s12 = s.col(j).head(j); // 获取s的第j列

        arma::vec beta = solveLasso(W11, s12, rho, s);

        Betas.col(j).head(j) = beta;  // 保存beta值到Betas矩阵

        arma::vec product = W11 * beta;  // 将乘积结果存储在向量中
        W.submat(0, j, j-1, j) = product;
        W.submat(j, 0, j, j-1) = product.t();
      }
        Rcpp::Rcout << "Completed iteration " << iter << " with W:\n" << W << "\n";

          // 检查收敛性
          if (maxDifference(W, W_old) < tol) {
            break;
          }

      arma::mat Theta = computeTheta(W, Betas);  // 调用函数计算\Theta，并传递Betas矩阵

      Rcpp::Rcout << "Completed glarma_cpp.\n";

      return List::create(Named("W") = W, Named("Theta") = Theta);  // 返回W和\Theta两个矩阵
  }
  }
