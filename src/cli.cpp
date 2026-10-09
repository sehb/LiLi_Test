#include "info_matrix.h"

#include <Eigen/Dense>
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

using Vec3 = Eigen::Vector3d;

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [--scene planar|random] [--n N] [--rel r] [--abs a]\n";
}

int main(int argc, char** argv) {
    std::string scene = "planar";
    int n = 200;
    double rel = 1e-3;
    double abs_th = 1e-9;

    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if (s == "--scene" && i + 1 < argc) { scene = argv[++i]; }
        else if (s == "--n" && i + 1 < argc) { n = std::atoi(argv[++i]); }
        else if (s == "--rel" && i + 1 < argc) { rel = std::atof(argv[++i]); }
        else if (s == "--abs" && i + 1 < argc) { abs_th = std::atof(argv[++i]); }
        else {
            print_usage(argv[0]);
            return 1;
        }
    }

    std::vector<std::tuple<Vec3, Vec3, Vec3>> data;
    if (scene == "planar") {
        for (int i = 0; i < n; ++i) {
            const double x = static_cast<double>(i) / n * 2.0 - 1.0;
            const double y = static_cast<double>((i * 37) % 101) / 50.0 - 1.0;
            const Vec3 p(x, y, 0.0);
            const Vec3 normal(0.0, 0.0, 1.0);
            data.emplace_back(p, p, normal);
        }
    } else if (scene == "random") {
        std::mt19937_64 rng(1234ULL);
        std::uniform_real_distribution<double> unif(-1.0, 1.0);
        std::uniform_real_distribution<double> ang(0.0, 2.0 * M_PI);
        for (int i = 0; i < n; ++i) {
            Vec3 p(unif(rng), unif(rng), unif(rng));
            const double theta = ang(rng);
            const double phi = ang(rng);
            Vec3 normal(std::cos(theta) * std::sin(phi),
                        std::sin(theta) * std::sin(phi),
                        std::cos(phi));
            normal.normalize();
            data.emplace_back(p, p, normal);
        }
    } else {
        std::cerr << "Unknown scene: " << scene << "\n";
        print_usage(argv[0]);
        return 1;
    }

    InformationMatrixCalculator calc;
    for (const auto& it : data) {
        Vec3 p, q, n;
        std::tie(p, q, n) = it;
        calc.addPointPlane(p, q, n);
    }

    const auto res = calc.analyzeDegeneracy(rel, abs_th);
    std::cout << "Scene: " << scene << " (N=" << n << ")\n";
    std::cout << "min_eigenvalue=" << res.min_eigenvalue << "\n";
    std::cout << "max_eigenvalue=" << res.max_eigenvalue << "\n";
    std::cout << "condition_number=" << res.condition_number << "\n";
    std::cout << "small_modes=" << res.small_indices.size() << "\n";
    for (size_t k = 0; k < res.small_indices.size(); ++k) {
        std::cout << "  mode=" << res.small_indices[k]
                  << " type=" << res.small_classification[k]
                  << " lambda=" << res.eigenvalues(res.small_indices[k]) << "\n";
    }
    return 0;
}
