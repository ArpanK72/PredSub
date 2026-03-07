#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;

// ---------------------------
  // Helper: build I_{p,q} via sign flips
// ---------------------------
  inline RowVectorXd build_sign_vector(int d, int p) {
    RowVectorXd s = RowVectorXd::Ones(d);
    if (p < d) {
      s.tail(d - p).setConstant(-1.0);
    }
    return s;
  }

// ---------------------------
  // 1) Full graph: sparsegraph_under_H0
//    R: sparsegraph_under_H0(X, Y, blocks, p)
// ---------------------------
  
  // [[Rcpp::export]]
SEXP sparsegraph_under_H0_rcpp(const MatrixXd& X,
                               const MatrixXd& Y,
                               int blocks,
                               int p) {
  const Eigen::Index n = X.rows();
  if (n != Y.rows() || X.cols() != Y.cols()) {
    Rcpp::stop("Matrices X and Y must have the same shape (n x d)");
  }
  const Eigen::Index d = X.cols();
  if (p < 0 || p > d) {
    Rcpp::stop("p must be in [0, d]");
  }
  
  // GRDPG sign pattern: implements I_{p,q}
  RowVectorXd s = build_sign_vector(d, p);
  
  // X*I_{p,q}, Y*I_{p,q} via rowwise sign flips
  MatrixXd XD = X.array().rowwise() * s.array();
  MatrixXd YD = Y.array().rowwise() * s.array();
  
  std::vector<Triplet<double> > trips;
  trips.reserve(static_cast<size_t>(n) * 20u); // heuristic similar to R code
  
  // Block over rows [1..n] as in your R code
  for (int blk = 1; blk <= blocks; ++blk) {
    
    // R-style block indices (1-based)
    double a_d = std::round((blk - 1) * (double)n / blocks + 1.0);
    double b_d = std::round( blk       * (double)n / blocks);
    int aR = static_cast<int>(a_d);
    int bR = static_cast<int>(b_d);
    if (bR > n) bR = n;
    if (aR > bR) continue;
    
    // Convert to 0-based for Eigen
    int a = aR - 1;
    int b = bR - 1;
    Eigen::Index n_rows = b - a + 1;
    
    // M = 0.5 * ( XD[a:b,] * X^T + YD[a:b,] * Y^T )  (size: (b-a+1) x n)
    MatrixXd XD_blk = XD.block(a, 0, n_rows, d);
    MatrixXd YD_blk = YD.block(a, 0, n_rows, d);
    MatrixXd M = 0.5 * (XD_blk * X.transpose() + YD_blk * Y.transpose());
    
    // For each i in [a,b], use row i-a of M and columns j>i
    for (int i = a; i <= b; ++i) {
      int num_trials = n - (i + 1);
      if (num_trials <= 0) continue;
      
      int rowM = i - a; // 0-based row index within M
      
      for (int j = i + 1; j < n; ++j) {
        double pij = M(rowM, j);
        // No explicit clamp, mimics your R code (p<0 => never, p>1 => always)
        if (R::unif_rand() < pij) {
          trips.emplace_back(i, j, 1.0);
          trips.emplace_back(j, i, 1.0);
        }
      }
    }
  }
  
  SparseMatrix<double> A(n, n);
  A.setFromTriplets(trips.begin(), trips.end());
  A.makeCompressed();
  
  return Rcpp::wrap(A); // dgCMatrix
}

// ---------------------------
  // 2) PredSub: sparsegraph_predsub_under_H0
//    R: sparsegraph_predsub_under_H0(X, Y, block, S, p)
// ---------------------------
  
  // [[Rcpp::export]]
Rcpp::List sparsegraph_predsub_under_H0_rcpp(const MatrixXd& X,
                                             const MatrixXd& Y,
                                             int block,
                                             Rcpp::IntegerVector S,
                                             int p) {
  const Eigen::Index n = X.rows();
  if (n != Y.rows() || X.cols() != Y.cols()) {
    Rcpp::stop("Matrices X and Y must have the same shape (n x d)");
  }
  const Eigen::Index d = X.cols();
  if (p < 0 || p > d) {
    Rcpp::stop("p must be in [0, d]");
  }
  
  // Build sign pattern and flipped matrices once
  RowVectorXd s = build_sign_vector(d, p);
  MatrixXd XD = X.array().rowwise() * s.array();
  MatrixXd YD = Y.array().rowwise() * s.array();
  
  // ---- Build S (0-based, unique, sorted) and R = complement ----
    std::vector<char> inS(static_cast<size_t>(n), 0);
  std::vector<Eigen::Index> S0;
  S0.reserve(S.size());
  
  for (int k = 0; k < S.size(); ++k) {
    int idxR = S[k];
    if (idxR < 1 || idxR > n) {
      Rcpp::stop("S must be a subset of {1,...,n}.");
    }
    Eigen::Index idx = idxR - 1; // 0-based
    if (!inS[idx]) {
      inS[idx] = 1;
      S0.push_back(idx);
    }
  }
  std::sort(S0.begin(), S0.end());
  
  std::vector<Eigen::Index> R0;
  R0.reserve(static_cast<size_t>(n) - S0.size());
  for (Eigen::Index i = 0; i < n; ++i) {
    if (!inS[i]) R0.push_back(i);
  }
  
  const Eigen::Index m         = static_cast<Eigen::Index>(S0.size());
  const Eigen::Index n_minus_m = static_cast<Eigen::Index>(R0.size());
  
  // ---------- Build B (S x S) using the full-graph function ----------
    // Gather X_s, Y_s (S-rows)
  MatrixXd Xs(m, d), Ys(m, d);
  for (Eigen::Index t = 0; t < m; ++t) {
    Xs.row(t) = X.row(S0[t]);
    Ys.row(t) = Y.row(S0[t]);
  }
  // Call the C++ full-graph function on the S-submatrix
  SEXP B_sexp = sparsegraph_under_H0_rcpp(Xs, Ys, block, p);
  Rcpp::S4 B(B_sexp); // dgCMatrix
  
  // ---------- Build A (R x S), rectangular ----------
    // We want probabilities:
    //   P_{ij} = 0.5 * (X_i^T I_{p,q} X_j + Y_i^T I_{p,q} Y_j)
  // for i in R, j in S. With XD = X I_{p,q}, YD = Y I_{p,q},
  // this is 0.5 * (XD_i^T X_j + YD_i^T Y_j).
  
  // Gather flipped rows for R and S
  MatrixXd XrD(n_minus_m, d), YrD(n_minus_m, d);
  MatrixXd Xs_orig(m, d), Ys_orig(m, d); // original X,Y rows at S
  for (Eigen::Index t = 0; t < n_minus_m; ++t) {
    XrD.row(t) = XD.row(R0[t]);
    YrD.row(t) = YD.row(R0[t]);
  }
  for (Eigen::Index t = 0; t < m; ++t) {
    Xs_orig.row(t) = X.row(S0[t]);
    Ys_orig.row(t) = Y.row(S0[t]);
  }
  
  std::vector<Triplet<double> > tripsA;
  tripsA.reserve(static_cast<size_t>(n_minus_m) * 10u); // heuristic
  
  // Loop over columns j_idx = 0..m-1 (S positions), as in your R code
  for (Eigen::Index j_idx = 0; j_idx < m; ++j_idx) {
    const RowVectorXd xj = Xs_orig.row(j_idx);
    const RowVectorXd yj = Ys_orig.row(j_idx);
    
    // p_vec = 0.5 * (X_r %*% xj + Y_r %*% yj)  where X_r = X[R,] %*% I_pq => XrD
    VectorXd p_vec = 0.5 * (XrD * xj.transpose() + YrD * yj.transpose());
    
    for (Eigen::Index i_r = 0; i_r < n_minus_m; ++i_r) {
      double pij = p_vec[i_r]; // no clamp, like your R
      if (R::unif_rand() < pij) {
        // A(i_r, j_idx) = 1
        tripsA.emplace_back(i_r, j_idx, 1.0);
      }
    }
  }
  
  SparseMatrix<double> A_sparse(n_minus_m, m);
  A_sparse.setFromTriplets(tripsA.begin(), tripsA.end());
  A_sparse.makeCompressed();
  
  Rcpp::S4 A(Rcpp::wrap(A_sparse)); // dgCMatrix
  
  return Rcpp::List::create(
    Rcpp::Named("B") = B,
    Rcpp::Named("A") = A
  );
}

#pragma GCC diagnostic pop
