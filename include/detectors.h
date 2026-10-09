#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "icp.h"
#include "point_cloud.h"
#include "se3.h"

namespace lili {

// A degeneracy subspace expressed as columns of se(3) twists in the *body*
// frame, matching the convention of paper Eq. (7) and (9):
//   T_i = P_opt * exp( sum_j c_ij * b_j ),   b_j = basis.col(j).
struct DetectorResult {
    std::string method;
    Eigen::Matrix<double, 6, Eigen::Dynamic> basis;
    Eigen::VectorXd eigenvalues;         // ascending, of the analyzed matrix
    Eigen::Matrix<double, 6, 6> eigenvectors;
    std::vector<int> degenerate_indices;
    double condition_number = 0.0;
    double threshold = 0.0;
};

// Zhang's Hessian-based degeneracy detection: eigen-decompose the information
// matrix and keep the directions whose eigenvalue falls below
// max(rel_threshold * max(1, lambda_max), abs_threshold).
struct ZhangOptions {
    double rel_threshold = 1e-3;
    double abs_threshold = 1e-9;
};

DetectorResult zhangDetector(const Eigen::Matrix<double, 6, 6>& information,
                             const ZhangOptions& options = ZhangOptions{});

}  // namespace lili
