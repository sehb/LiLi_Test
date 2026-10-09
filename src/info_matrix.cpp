#include "info_matrix.h"

#include <Eigen/Eigenvalues>
#include <cmath>

namespace {
Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
    Eigen::Matrix3d m;
    m << 0.0, -v.z(), v.y(),
         v.z(), 0.0, -v.x(),
         -v.y(), v.x(), 0.0;
    return m;
}
}

InformationMatrixCalculator::InformationMatrixCalculator() {
    reset();
}

void InformationMatrixCalculator::reset() {
    info_.setZero();
}

void InformationMatrixCalculator::addPointPlane(const Eigen::Vector3d& p_src,
                                               const Eigen::Vector3d& q_tgt,
                                               const Eigen::Vector3d& n_tgt,
                                               const Eigen::Matrix3d& R,
                                               const Eigen::Vector3d& t) {
    (void)q_tgt;
    (void)t;

    const Eigen::RowVector3d dr_dphi = -(n_tgt.transpose() * R * skew(p_src));
    const Eigen::RowVector3d dr_drho = n_tgt.transpose();

    Eigen::Matrix<double, 1, 6> J;
    J << dr_dphi, dr_drho;

    info_ += J.transpose() * J;
}

Eigen::Matrix<double, 6, 6> InformationMatrixCalculator::informationMatrix() const {
    return info_;
}

DegeneracyResult InformationMatrixCalculator::analyzeDegeneracy(double rel_threshold,
                                                                double abs_threshold) const {
    DegeneracyResult res;
    res.information = info_;

    Eigen::Matrix<double, 6, 6> sym = (info_ + info_.transpose()) * 0.5;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> solver(sym);

    res.eigenvalues = solver.eigenvalues();
    res.eigenvectors = solver.eigenvectors();

    const double lambda_max = res.eigenvalues.maxCoeff();
    const double rel_cut = rel_threshold * std::max(1.0, lambda_max);

    for (int i = 0; i < 6; ++i) {
        const double lambda = res.eigenvalues(i);
        if (lambda < rel_cut || lambda < abs_threshold) {
            res.small_indices.push_back(i);

            const Eigen::VectorXd v = res.eigenvectors.col(i);
            const double rot_norm = v.head<3>().norm();
            const double trans_norm = v.tail<3>().norm();

            std::string cls;
            if (rot_norm > 2.0 * trans_norm) {
                cls = "rotation";
            } else if (trans_norm > 2.0 * rot_norm) {
                cls = "translation";
            } else {
                cls = "mixed";
            }
            res.small_classification.push_back(cls);
        }
    }

    return res;
}
