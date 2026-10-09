#include "detectors.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace lili {

DetectorResult zhangDetector(const Eigen::Matrix<double, 6, 6>& information,
                             const ZhangOptions& options) {
    DetectorResult result;
    result.method = "zhang_hessian";

    const Eigen::Matrix<double, 6, 6> sym = 0.5 * (information + information.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> solver(sym);
    result.eigenvalues = solver.eigenvalues();      // ascending
    result.eigenvectors = solver.eigenvectors();

    const double lambda_max = result.eigenvalues.maxCoeff();
    const double lambda_min = result.eigenvalues.minCoeff();
    result.condition_number =
        lambda_min > 0.0 ? lambda_max / lambda_min
                         : std::numeric_limits<double>::infinity();

    const double cut = std::max(options.rel_threshold * std::max(1.0, lambda_max),
                                options.abs_threshold);
    result.threshold = cut;

    std::vector<int> indices;
    for (int i = 0; i < 6; ++i) {
        if (result.eigenvalues(i) < cut) { indices.push_back(i); }
    }

    result.basis.resize(6, static_cast<int>(indices.size()));
    for (std::size_t k = 0; k < indices.size(); ++k) {
        result.basis.col(static_cast<int>(k)) = result.eigenvectors.col(indices[k]);
    }
    result.degenerate_indices = indices;
    return result;
}

}  // namespace lili
