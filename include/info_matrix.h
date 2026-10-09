#pragma once

#include <Eigen/Dense>
#include <limits>
#include <string>
#include <vector>

struct DegeneracyResult {
    Eigen::Matrix<double, 6, 6> information;
    Eigen::VectorXd eigenvalues;
    Eigen::Matrix<double, 6, 6> eigenvectors;
    std::vector<int> small_indices;
    std::vector<std::string> small_classification;
    double min_eigenvalue = 0.0;
    double max_eigenvalue = 0.0;
    double condition_number = std::numeric_limits<double>::infinity();
};

class InformationMatrixCalculator {
public:
    InformationMatrixCalculator();

    void reset();

    // Add one point-to-plane residual constraint:
    // residual r = n^T * (R * p_src + t - q_tgt)
    // with Jacobian
    // J = [ n^T * (-R * [p_src]_x) , n^T ]
    void addPointPlane(const Eigen::Vector3d& p_src,
                       const Eigen::Vector3d& q_tgt,
                       const Eigen::Vector3d& n_tgt,
                       const Eigen::Matrix3d& R = Eigen::Matrix3d::Identity(),
                       const Eigen::Vector3d& t = Eigen::Vector3d::Zero());

    Eigen::Matrix<double, 6, 6> informationMatrix() const;

    // Lower eigenvalues indicate a weakly observable direction.
    DegeneracyResult analyzeDegeneracy(double rel_threshold = 1e-3,
                                      double abs_threshold = 1e-9) const;

    // Condition number of the accumulated information matrix.
    double conditionNumber() const;

private:
    Eigen::Matrix<double, 6, 6> info_;
};
