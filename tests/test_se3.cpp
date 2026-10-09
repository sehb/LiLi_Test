#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>

#include "se3.h"

namespace {

int failures = 0;

void expect(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    }
}

}  // namespace

int main() {
    using namespace lili;

    // hat/vee round trip and nilpotency of the twist matrix.
    {
        Vector6d xi;
        xi << 0.3, -0.2, 0.5, 1.0, -2.0, 0.25;
        expect((vee(hat(xi)) - xi).norm() < 1e-12, "hat/vee round trip");
        const Eigen::Matrix4d X = hat(xi);
        expect(X.bottomRightCorner<1, 1>().norm() < 1e-12,
               "hat of a twist is nilpotent");
    }

    // exp/log round trip. The log is unique only for |phi| < pi, so the
    // rotation magnitude is kept strictly inside that range.
    {
        std::mt19937_64 rng(7u);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        double worst = 0.0;
        for (int trial = 0; trial < 200; ++trial) {
            Eigen::Vector3d axis(u(rng), u(rng), u(rng));
            if (axis.norm() < 1e-6) { axis = Eigen::Vector3d::UnitX(); }
            Vector6d xi;
            xi.head<3>() = axis.normalized() * (0.02 + (M_PI - 0.06) * std::abs(u(rng)));
            xi.tail<3>() = Eigen::Vector3d(u(rng), u(rng), u(rng));
            const SE3 T = expSE3(xi);
            worst = std::max(worst, (logSE3(T) - xi).norm());
        }
        expect(worst < 1e-8, "exp/log round trip");
    }

    // For large twists log returns the principal value, so compare through the
    // group instead: exp(log(T)) must reproduce T exactly.
    {
        std::mt19937_64 rng(11u);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        double worst = 0.0;
        for (int trial = 0; trial < 200; ++trial) {
            Vector6d xi;
            xi.head<3>() = Eigen::Vector3d(u(rng), u(rng), u(rng)) * 3.0 * M_PI;
            xi.tail<3>() = Eigen::Vector3d(u(rng), u(rng), u(rng)) * 2.0;
            const SE3 T = expSE3(xi);
            const SE3 round = expSE3(logSE3(T));
            worst = std::max(worst, (T.R - round.R).norm() + (T.t - round.t).norm());
        }
        expect(worst < 1e-9, "exp(log(T)) == T for large twists");
    }

    // exp produces a valid rotation.
    {
        Vector6d xi;
        xi << 0.1, 0.2, -0.3, 0.5, 0.0, -0.7;
        const SE3 T = expSE3(xi);
        expect((T.R.transpose() * T.R - Eigen::Matrix3d::Identity()).norm() < 1e-12,
               "exp produces an orthonormal R");
        expect(std::abs(T.R.determinant() - 1.0) < 1e-12, "det(R) == 1");
    }

    // Composition and inverse.
    {
        Vector6d a, b;
        a << 0.2, -0.1, 0.3, 1.0, 0.5, -0.2;
        b << -0.4, 0.25, 0.1, -0.3, 0.8, 0.6;
        const SE3 Ta = expSE3(a), Tb = expSE3(b);
        const Eigen::Vector3d p(0.4, -0.2, 0.9);
        expect(((Ta * Tb) * p - Ta * (Tb * p)).norm() < 1e-12,
               "composition acts as R p + t");
        const SE3 T = Ta * Tb;
        expect((T * T.inverse()).t.norm() < 1e-12, "inverse composition");
        expect((T.R * T.R.transpose() - Eigen::Matrix3d::Identity()).norm() < 1e-12,
               "R stays orthogonal");
    }

    // Adjoint: Ad_T Ad_{T^-1} == I and conjugation T exp(xi) T^-1 == exp(Ad_T xi).
    {
        Vector6d a, xi;
        a << 0.3, 0.2, -0.1, 0.7, -0.4, 0.2;
        xi << -0.2, 0.5, 0.15, 0.3, 0.1, -0.6;
        const SE3 T = expSE3(a);
        const Matrix6d Ad = adjoint(T);
        expect((Ad * adjoint(T.inverse()) - Matrix6d::Identity()).norm() < 1e-10,
               "Ad_T Ad_{T^-1} == I");
        const SE3 lhs = T * expSE3(xi) * T.inverse();
        const SE3 rhs = expSE3(Ad * xi);
        expect((lhs.R - rhs.R).norm() < 1e-9, "conjugation rotation");
        expect((lhs.t - rhs.t).norm() < 1e-9, "conjugation translation");
    }

    // twistAtPoint moves the linear part by phi x offset.
    {
        Vector6d xi = Vector6d::Zero();
        xi.head<3>() = Eigen::Vector3d::UnitZ();
        const Eigen::Vector3d sensor(1.0, 0.0, 0.0);
        const Vector6d moved = twistAtPoint(xi, sensor);
        expect((moved.head<3>() - xi.head<3>()).norm() < 1e-12,
               "twistAtPoint keeps the angular part");
        expect((moved.tail<3>() - Eigen::Vector3d(0.0, 1.0, 0.0)).norm() < 1e-12,
               "twistAtPoint offsets the linear part by phi x c");
    }

    if (failures == 0) {
        std::cout << "LiLi SE(3) tests passed.\n";
    }
    return failures == 0 ? 0 : 1;
}
