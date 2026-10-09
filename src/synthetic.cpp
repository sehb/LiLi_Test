#include "info_matrix.h"

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <tuple>
#include <vector>

using Vec3 = Eigen::Vector3d;

static std::mt19937_64 rng(123456789ULL);

std::vector<std::tuple<Vec3, Vec3, Vec3>> makePlanarScene(int n) {
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::vector<std::tuple<Vec3, Vec3, Vec3>> out;
    out.reserve(n);

    const Vec3 normal(0.0, 0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        const double x = unif(rng);
        const double y = unif(rng);
        const Vec3 p_src(x, y, 0.0);
        const Vec3 q_tgt = p_src;
        out.emplace_back(p_src, q_tgt, normal);
    }
    return out;
}

std::vector<std::tuple<Vec3, Vec3, Vec3>> makeRandom3DScene(int n) {
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::uniform_real_distribution<double> ang(0.0, 2.0 * M_PI);
    std::vector<std::tuple<Vec3, Vec3, Vec3>> out;
    out.reserve(n);

    for (int i = 0; i < n; ++i) {
        const Vec3 p_src(unif(rng), unif(rng), unif(rng));
        const Vec3 q_tgt = p_src;
        const double theta = ang(rng);
        const double phi = ang(rng);
        Vec3 normal(std::cos(theta) * std::sin(phi),
                    std::sin(theta) * std::sin(phi),
                    std::cos(phi));
        normal.normalize();
        out.emplace_back(p_src, q_tgt, normal);
    }
    return out;
}

void printDegeneracy(const std::string& name, const DegeneracyResult& r) {
    std::cout << "Scene: " << name << "\n";
    std::cout << "Eigenvalues (ascending):";
    for (int i = 0; i < 6; ++i) {
        std::cout << " " << std::setw(12) << r.eigenvalues(i);
    }
    std::cout << "\n";

    if (r.small_indices.empty()) {
        std::cout << "No degenerate modes detected.\n";
    } else {
        std::cout << "Detected small eigenvalues:\n";
        for (size_t k = 0; k < r.small_indices.size(); ++k) {
            const int idx = r.small_indices[k];
            std::cout << " - index=" << idx
                      << ", lambda=" << r.eigenvalues(idx)
                      << ", type=" << r.small_classification[k] << "\n";
        }
    }
    std::cout << "--------------------------------\n";
}

int main() {
    {
        auto data = makePlanarScene(200);
        InformationMatrixCalculator calc;
        for (const auto& item : data) {
            Vec3 p, q, n;
            std::tie(p, q, n) = item;
            calc.addPointPlane(p, q, n);
        }
        const auto result = calc.analyzeDegeneracy(1e-3, 1e-9);
        printDegeneracy("Planar", result);
        if (result.small_indices.empty()) {
            std::cerr << "Planar scene should be degenerate!\n";
            return 1;
        }
    }

    {
        auto data = makeRandom3DScene(500);
        InformationMatrixCalculator calc;
        for (const auto& item : data) {
            Vec3 p, q, n;
            std::tie(p, q, n) = item;
            calc.addPointPlane(p, q, n);
        }
        const auto result = calc.analyzeDegeneracy(1e-3, 1e-9);
        printDegeneracy("Random3D", result);
        if (!result.small_indices.empty()) {
            std::cerr << "Random 3D scene unexpectedly degenerate.\n";
            return 1;
        }
    }

    std::cout << "Synthetic validation passed.\n";
    return 0;
}
