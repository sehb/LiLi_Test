#include "se3.h"

#include <cmath>

namespace lili {
namespace {

// Extract [v]_x -> v from the skew-symmetric part of a 3x3 matrix.
Eigen::Vector3d vee3(const Eigen::Matrix3d& M) {
    return Eigen::Vector3d(M(2, 1), M(0, 2), M(1, 0));
}

// V(phi) = I + ((1-cos t)/t^2) [phi]_x + ((t-sin t)/t^3) [phi]_x^2
Eigen::Matrix3d leftJacobian(const Eigen::Vector3d& phi) {
    const double t2 = phi.squaredNorm();
    const Eigen::Matrix3d phix = skew(phi);
    if (t2 < 1e-8) {
        return Eigen::Matrix3d::Identity() + 0.5 * phix + (1.0 / 6.0) * phix * phix;
    }
    const double t = std::sqrt(t2);
    return Eigen::Matrix3d::Identity() + ((1.0 - std::cos(t)) / t2) * phix +
           ((t - std::sin(t)) / (t2 * t)) * phix * phix;
}

// V(phi)^-1 = I - 0.5 [phi]_x + (1/t^2)(1 - t sin t / (2(1-cos t))) [phi]_x^2
Eigen::Matrix3d leftJacobianInverse(const Eigen::Vector3d& phi) {
    const double t2 = phi.squaredNorm();
    const Eigen::Matrix3d phix = skew(phi);
    if (t2 < 1e-8) {
        return Eigen::Matrix3d::Identity() - 0.5 * phix + (1.0 / 12.0) * phix * phix;
    }
    const double t = std::sqrt(t2);
    const double coeff =
        (1.0 - (t * std::sin(t)) / (2.0 * (1.0 - std::cos(t)))) / t2;
    return Eigen::Matrix3d::Identity() - 0.5 * phix + coeff * phix * phix;
}

}  // namespace

SE3 SE3::identity() { return SE3{}; }

SE3 SE3::inverse() const {
    SE3 inv;
    inv.R = R.transpose();
    inv.t = -inv.R * t;
    return inv;
}

Eigen::Vector3d SE3::operator*(const Eigen::Vector3d& p) const {
    return R * p + t;
}

SE3 SE3::operator*(const SE3& other) const {
    SE3 out;
    out.R = R * other.R;
    out.t = R * other.t + t;
    return out;
}

Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
    Eigen::Matrix3d m;
    m << 0.0, -v.z(), v.y(),
         v.z(), 0.0, -v.x(),
         -v.y(), v.x(), 0.0;
    return m;
}

Eigen::Matrix4d hat(const Vector6d& xi) {
    Eigen::Matrix4d X = Eigen::Matrix4d::Zero();
    X.topLeftCorner<3, 3>() = skew(xi.head<3>());
    X.topRightCorner<3, 1>() = xi.tail<3>();
    return X;
}

Vector6d vee(const Eigen::Matrix4d& X) {
    Vector6d xi;
    xi.head<3>() = vee3(X.topLeftCorner<3, 3>());
    xi.tail<3>() = X.topRightCorner<3, 1>();
    return xi;
}

SE3 expSE3(const Vector6d& xi) {
    const Eigen::Vector3d phi = xi.head<3>();
    const Eigen::Vector3d rho = xi.tail<3>();

    SE3 T;
    const double t2 = phi.squaredNorm();
    const Eigen::Matrix3d phix = skew(phi);
    if (t2 < 1e-8) {  // Rodrigues with the small-angle series
        T.R = Eigen::Matrix3d::Identity() + phix + 0.5 * phix * phix;
    } else {
        const double t = std::sqrt(t2);
        T.R = Eigen::Matrix3d::Identity() + (std::sin(t) / t) * phix +
              ((1.0 - std::cos(t)) / t2) * phix * phix;
    }
    T.t = leftJacobian(phi) * rho;
    return T;
}

Vector6d logSE3(const SE3& T) {
    const double c = std::max(-1.0, std::min(1.0, (T.R.trace() - 1.0) * 0.5));
    const double theta = std::acos(c);

    Eigen::Vector3d phi = Eigen::Vector3d::Zero();
    if (theta < 1e-10) {
        phi = 0.5 * vee3(T.R - T.R.transpose());
    } else if (theta < M_PI - 1e-6) {
        phi = (theta / (2.0 * std::sin(theta))) * vee3(T.R - T.R.transpose());
    } else {
        const Eigen::Matrix3d A = 0.5 * (T.R + Eigen::Matrix3d::Identity());
        Eigen::Vector3d axis = A.col(0);
        int best = 0;
        for (int i = 1; i < 3; ++i) {
            if (A(i, i) > A(best, best)) { best = i; }
        }
        axis = A.col(best);
        phi = theta * axis.normalized();
    }

    Vector6d xi;
    xi.head<3>() = phi;
    xi.tail<3>() = leftJacobianInverse(phi) * T.t;
    return xi;
}

Matrix6d adjoint(const SE3& T) {
    Matrix6d Ad = Matrix6d::Zero();
    Ad.topLeftCorner<3, 3>() = T.R;
    Ad.topRightCorner<3, 3>() = Eigen::Matrix3d::Zero();
    Ad.bottomLeftCorner<3, 3>() = skew(T.t) * T.R;
    Ad.bottomRightCorner<3, 3>() = T.R;
    return Ad;
}

Vector6d twistAtPoint(const Vector6d& xi, const Eigen::Vector3d& to,
                      const Eigen::Vector3d& from) {
    Vector6d out = xi;
    out.tail<3>() = xi.tail<3>() + xi.head<3>().cross(to - from);
    return out;
}

}  // namespace lili
