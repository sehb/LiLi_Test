#pragma once

#include <Eigen/Dense>
#include <vector>

#include "point_cloud.h"
#include "se3.h"

namespace lili {

// Extended-scan quality metric from the paper (Eq. 8-10):
//   T_i = P_opt * exp(sum_j c_ij * b_j),  S_ext = union_i S_local * T_i
//   Q   = median over points of S_ext of the NN distance to S_reference.
//
// The coefficients are drawn from a deterministic low-discrepancy sequence
// (Halton) so that two methods being compared receive *identical* coefficients
// even when their detected bases have different dimensions.
struct ExtensionOptions {
    int n_samples = 32;         // number of degenerate motions sampled
    double coefficient_scale = 1.0;  // max |c_ij|
    // Nearest-neighbour search cap. Q is a distance to the reference surface,
    // so it only has to discriminate "on the surface" (~one spacing) from "off
    // it". Capping makes Q bounded and keeps the grid search cheap; < 0 selects
    // 8 x the reference's median nearest-neighbour spacing.
    double max_radius = -1.0;
};

struct QualityResult {
    double median = 0.0;
    double mean = 0.0;
    double p90 = 0.0;
    int samples = 0;
};

QualityResult alignmentQuality(
    const PointCloud& local, const PointCloud& reference, const SE3& P_opt,
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& basis,
    const ExtensionOptions& opt = ExtensionOptions{});

// Convenience: Q of a single transform, i.e. no degeneracy extension.
QualityResult alignmentQualitySingle(const PointCloud& local,
                                     const PointCloud& reference,
                                     const SE3& T);

// Principal angles (degrees, descending) between the column spaces of two
// se(3) bases. Rotation entries are scaled by `length_scale` first so that the
// mixed rad/metre units of a twist do not distort the comparison.
std::vector<double> principalAnglesDeg(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& A,
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& B, double length_scale = 1.0);

// Decomposition of a body-frame twist against a degeneracy basis. Rotation
// entries are scaled by `length_scale` first, so the split is unit-consistent.
struct SubspaceDeviation {
    double total = 0.0;       // ||twist||
    double orthogonal = 0.0;  // ||twist - projection onto span(basis)||
    double fraction = 0.0;    // orthogonal / total (0 when total == 0)
};

SubspaceDeviation deviationFromSubspace(
    const Vector6d& twist,
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& basis,
    double length_scale = 1.0);

}  // namespace lili
