#pragma once

#include <Eigen/Dense>

namespace lili {

using Vector6d = Eigen::Matrix<double, 6, 1>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;

// Rigid transform, local -> reference: p_ref = R * p_local + t.
struct SE3 {
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    Eigen::Vector3d t = Eigen::Vector3d::Zero();

    static SE3 identity();
    SE3 inverse() const;
    Eigen::Vector3d operator*(const Eigen::Vector3d& p) const;
    SE3 operator*(const SE3& other) const;
};

Eigen::Matrix3d skew(const Eigen::Vector3d& v);

// Twist ordering matches info_matrix.h: xi = [phi (rotation); rho (translation)].
Eigen::Matrix4d hat(const Vector6d& xi);
Vector6d vee(const Eigen::Matrix4d& X);

SE3 expSE3(const Vector6d& xi);
Vector6d logSE3(const SE3& T);

// Ad_T maps a twist expressed in the body frame to the reference frame.
Matrix6d adjoint(const SE3& T);

// Re-express a twist given at `from` for a frame whose origin sits at `to`,
// both frames sharing the same orientation: rho_to = rho_from + phi x (to - from).
Vector6d twistAtPoint(const Vector6d& xi, const Eigen::Vector3d& to,
                      const Eigen::Vector3d& from = Eigen::Vector3d::Zero());

}  // namespace lili
