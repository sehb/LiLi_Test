#include "info_matrix.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

int main() {
    InformationMatrixCalculator calc;
    for (int i = 0; i < 200; ++i) {
        const double x = static_cast<double>(i) / 200.0;
        const double y = std::sin(static_cast<double>(i) * 0.37);
        const double z = std::cos(static_cast<double>(i) * 0.91);

        Eigen::Vector3d p(x, y, z);
        Eigen::Vector3d q = p;
        Eigen::Vector3d n(1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0));
        calc.addPointPlane(p, q, n);
    }

    const auto result = calc.analyzeDegeneracy(1e-3, 1e-9);
    const auto info = result.information;
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> solver(info);
    const auto lambdas = solver.eigenvalues();

    if (lambdas.minCoeff() < 0.0) {
        std::cerr << "Information matrix not positive semidefinite.\n";
        return 1;
    }

    std::cout << "LiLi info-matrix tests passed.\n";
    return 0;
}
