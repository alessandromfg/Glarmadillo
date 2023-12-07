// glarma.cpp
// Part of the Glarmadillo package
// Copyright (C) 2023 [Alessandro]
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// 软阈值计算
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
arma::vec solveLasso(const arma::mat& V, const arma::vec& u, double rho, const arma::mat& s, double ltol) {
  int p = V.n_rows;
  arma::vec beta = arma::zeros(p);
  arma::vec beta_old(p);
  arma::vec off_diag;

  off_diag = s.elem(arma::find(arma::trimatu(arma::ones<arma::mat>(p,p), 1)));  // 获取非主对角线的元素

  double threshold = ltol * arma::mean(arma::abs(off_diag));  // 计算阈值
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


double maxDifference(const arma::mat& A, const arma::mat& B) {
  return arma::max(arma::max(arma::abs(A - B)));
}

// 计算Theta矩阵的函数
arma::mat computeTheta(const arma::mat& W, const arma::mat& Betas) {
  int p = W.n_rows;
  arma::mat Theta(p, p, arma::fill::zeros); // 初始化Theta为p x p的零矩阵

  for (int j = 0; j < p; j++) {
    arma::vec w12 = W.col(j);
    w12.shed_row(j); // 去掉对角线元素的第j列
    double w22 = W(j, j);
    arma::vec beta = Betas.col(j); // beta是一个p-1维的向量

    double theta22 = 1.0 / (w22 - arma::dot(w12, beta));
    arma::vec theta12 = -beta * theta22; // theta12是一个p-1维的向量

    // 将theta22赋值到Theta矩阵的对角线位置
    Theta(j, j) = theta22;

    if (j == 0) {
      // 当j为0时，只更新Theta矩阵的下半部分
      Theta.submat(j + 1, j, p - 1, j) = theta12;
    } else if (j == p - 1) {
      // 当j为p-1时，只更新Theta矩阵的上半部分
      Theta.submat(0, j, j - 1, j) = theta12;
    } else {
      // 其他情况，更新Theta矩阵的上半部分和下半部分
      Theta.submat(0, j, j - 1, j) = theta12.head(j);
      Theta.submat(j + 1, j, p - 1, j) = theta12.subvec(j, p - 2);
    }
  }
  return Theta; // 返回计算得到的\Theta矩阵
}

// 用于 n=1 的情况的特殊处理函数
List solveGlarmaForOne(const arma::mat& s, double rho) {
  double W = s(0,0) + rho;  // 因为 s 是一个标量，这里 W 也是一个标量
  double Theta = 1.0 / W;   // 同样，Theta 也是一个标量

  // 构造输出结果
  return List::create(
    Named("W") = W,
    Named("Theta") = Theta,
    Named("converged") = true,
    Named("maxDiff") = 0.0,
    Named("iter") = 1  // 这里没有迭代，但为了一致性，设置为 1
  );
}

// 专门用于处理 n=2 的情况
List solveGlarmaForTwo(const arma::mat& s, double rho, double ltol) {
  // 直接计算 W 矩阵
  arma::mat W = s + rho * arma::eye(2, 2);

  // 计算 beta，只有一个元素
  double beta = softThreshold(s(0, 1), rho) / W(0, 0);

  // 检查 beta 的收敛性
  double beta_old;
  double threshold = ltol * std::abs(s(0, 1));
  do {
    beta_old = beta;
    beta = softThreshold(s(0, 1), rho) / W(0, 0);
  } while (std::abs(beta - beta_old) > threshold);

  // 更新 W 矩阵的非对角线元素
  W(0, 1) = W(0, 0) * beta;
  W(1, 0) = W(0, 1); // 对称矩阵

  // 计算 Theta 矩阵
  arma::mat Theta;
  try {
    Theta = W.i(); // 直接求 W 的逆矩阵
  } catch (std::runtime_error &e) {
    Rcpp::Rcout << "Error inverting W: " << e.what() << std::endl;
    Theta = arma::mat(2, 2, arma::fill::zeros); // 如果求逆失败，返回零矩阵
  }

  // 构造输出结果
  return List::create(
    Named("W") = W,
    Named("Theta") = Theta,
    Named("converged") = true,
    Named("maxDiff") = std::abs(beta - beta_old),
    Named("iter") = 1  // 这里没有迭代，但为了一致性，设置为 1
  );
}



// 主函数
// [[Rcpp::export]]
List glarma_cpp(const arma::mat& s, double rho, double mtol, int maxIterations, double ltol) {
  int p = s.n_rows;
  // 检查矩阵大小，对于n=1的情况进行特殊处理
  if (p == 1) {
    return solveGlarmaForOne(s, rho);
  }
  // 检查矩阵大小，对于n=2的情况进行特殊处理
  if (p == 2) {
    return solveGlarmaForTwo(s, rho, ltol);
  }

  arma::mat W = initializeW(s, rho);
  arma::mat W_old;

  bool converged = false; // 是否收敛的标志
  double maxDiff = 0; // 最大差异
  int iter; // 实际迭代次数

  // 初始化Betas矩阵为p-1 x p的零矩阵
  arma::mat Betas(p-1, p, arma::fill::zeros);

  for (iter = 0; iter < maxIterations; iter++) {

    W_old = W; // 保存W矩阵

    for (int j = 0; j < p; j++) {

      arma::mat W11 = W;
      W11.shed_row(j);
      W11.shed_col(j);

      arma::vec s12 = s.col(j);
      s12.shed_row(j);

      arma::vec beta = solveLasso(W11, s12, rho, s, ltol);

      Betas.col(j) = beta;

      arma::vec product = W11 * beta;

      // 更新W矩阵的对应部分
      if (j == 0) {
        W.submat(1, 0, p-1, 0) = product;// 更新第0列，从第1行到最后
        W.submat(0, 1, 0, p-1) = product.t();// 更新第0行，从第1列到最后
      } else if (j == p-1) {
        W.submat(0, p-1, p-2, p-1) = product;// 更新最后一列，从第0行到第p-2行
        W.submat(p-1, 0, p-1, p-2) = product.t();// 更新最后一行，从第0列到第p-2列
      } else {
        W.submat(0, j, j-1, j) = product.head(j);  // 更新W的第j列，直到第j-1行
        W.submat(j+1, j, p-1, j) = product.subvec(j, p-2);  // 更新W的第j列，从第j+1行到最后
        W.submat(j, 0, j, j-1) = product.head(j).t();  // 更新W的第j行，直到第j-1列
        W.submat(j, j+1, j, p-1) = product.subvec(j, p-2).t();  // 更新W的第j行，从第j+1列到最后
      }

    }

    // 检查收敛性
    maxDiff = maxDifference(W, W_old);
    if (maxDiff < mtol) {
      converged = true;
      break;
    }
  }

  arma::mat Theta = computeTheta(W, Betas); // 计算Theta
  // 构造输出结果
  return List::create(
    Named("W") = W,
    Named("Theta") = Theta,
    Named("converged") = converged,
    Named("maxDiff") = maxDiff,
    Named("iter") = iter + 1 // 实际迭代次数，从1开始计数
  );
}
