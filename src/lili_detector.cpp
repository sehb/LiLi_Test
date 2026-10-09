#include "lili_detector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace lili {
namespace {

// Median point displacement ||exp(xi) p - p|| produced by the body-frame twist
// xi applied to every point of `cloud`.
double medianTwistDisplacement(const Vector6d& xi, const PointCloud& cloud) {
    if (cloud.empty()) return 0.0;
    const SE3 T = expSE3(xi);
    std::vector<double> d;
    d.reserve(cloud.size());
    for (const auto& p : cloud.points) {
        d.push_back((T * p - p).norm());
    }
    std::sort(d.begin(), d.end());
    return d[d.size() / 2];
}

// Median distance of the cloud's points to their centroid, the "benchmark
// radius" the paper calibrates the perturbation against.
double medianRadius(const PointCloud& cloud) {
    if (cloud.empty()) return 0.0;
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    for (const auto& p : cloud.points) c += p;
    c /= static_cast<double>(cloud.size());
    std::vector<double> d;
    d.reserve(cloud.size());
    for (const auto& p : cloud.points) d.push_back((p - c).norm());
    std::sort(d.begin(), d.end());
    return d[d.size() / 2];
}

Eigen::Matrix<double, 6, Eigen::Dynamic> scaleBasis(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& B, double length_scale) {
    Eigen::Matrix<double, 6, Eigen::Dynamic> S = B;
    if (length_scale != 0.0) S.topRows<3>() *= length_scale;
    return S;
}

Eigen::Matrix<double, 6, Eigen::Dynamic> unscaleBasis(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& S, double length_scale) {
    Eigen::Matrix<double, 6, Eigen::Dynamic> B = S;
    if (length_scale != 0.0) B.topRows<3>() /= length_scale;
    return B;
}

// Given two columns (i, j) of B, find the Givens angle in [0, pi) that
// minimises the combined L1 norm of the rotated pair. A dense scan seeds a
// ternary-search refinement; this is the 2D subproblem of Eq. (14).
double bestGivensAngle(const Eigen::MatrixXd& B, int i, int j) {
    const Eigen::VectorXd a = B.col(i);
    const Eigen::VectorXd b = B.col(j);
    auto costAt = [&](double theta) {
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        return (c * a - s * b).cwiseAbs().sum() + (s * a + c * b).cwiseAbs().sum();
    };

    const int n_samples = 180;
    double best_theta = 0.0;
    double best_cost = std::numeric_limits<double>::infinity();
    for (int s = 0; s < n_samples; ++s) {
        const double theta = M_PI * s / n_samples;
        const double cost = costAt(theta);
        if (cost < best_cost) {
            best_cost = cost;
            best_theta = theta;
        }
    }

    double lo = best_theta - M_PI / n_samples;
    double hi = best_theta + M_PI / n_samples;
    for (int it = 0; it < 60; ++it) {
        const double m1 = lo + (hi - lo) / 3.0;
        const double m2 = hi - (hi - lo) / 3.0;
        if (costAt(m1) < costAt(m2)) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    return 0.5 * (lo + hi);
}

}  // namespace

PerturbationScale calibratePerturbation(const PointCloud& local,
                                        const PointCloud& reference,
                                        int k) {
    PerturbationScale scale;
    scale.spacing = medianNearestNeighborSpacing(reference, 400);
    scale.benchmark_radius = medianRadius(local);
    const int kk = std::max(1, k);
    // A translation of k spacings slides the benchmark point k samples; a
    // rotation moving a point at the benchmark radius by the same arc gives the
    // matching angular magnitude.
    scale.translation = kk * scale.spacing;
    scale.rotation = (scale.benchmark_radius > 1e-9)
                         ? kk * scale.spacing / scale.benchmark_radius
                         : kk * scale.spacing;
    return scale;
}

std::vector<Vector6d> generatePerturbations(const PerturbationScale& scale,
                                            bool bidirectional) {
    std::vector<Vector6d> dirs;
    for (int axis = 0; axis < 3; ++axis) {
        Vector6d v = Vector6d::Zero();
        v(axis) = scale.rotation;         // rotation about the axis (rad)
        v(3 + axis) = scale.translation;  // translation along the axis (m)
        dirs.push_back(v);
        if (bidirectional) dirs.push_back(-v);
    }
    return dirs;
}

Eigen::Matrix<double, 6, Eigen::Dynamic> l1SparsifyBasis(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& B_in,
    double length_scale,
    int max_iter) {
    if (B_in.cols() <= 0) return B_in;
    const int k = static_cast<int>(B_in.cols());
    const Eigen::MatrixXd scaled = scaleBasis(B_in, length_scale);
    Eigen::MatrixXd B =
        Eigen::MatrixXd(scaled.householderQr().householderQ()).leftCols(k);

    for (int iter = 0; iter < max_iter; ++iter) {
        bool changed = false;
        for (int i = 0; i < k; ++i) {
            for (int j = i + 1; j < k; ++j) {
                const double cost0 =
                    B.col(i).cwiseAbs().sum() + B.col(j).cwiseAbs().sum();
                const double theta = bestGivensAngle(B, i, j);
                const double c = std::cos(theta);
                const double s = std::sin(theta);
                const Eigen::VectorXd a = B.col(i);
                const Eigen::VectorXd b = B.col(j);
                const Eigen::VectorXd na = c * a - s * b;
                const Eigen::VectorXd nb = s * a + c * b;
                const double cost1 =
                    na.cwiseAbs().sum() + nb.cwiseAbs().sum();
                if (cost0 - cost1 > 1e-12) {
                    B.col(i) = na;
                    B.col(j) = nb;
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }

    return unscaleBasis(B, length_scale);
}

DetectorResult liliDetector(const PointCloud& local,
                            const PointCloud& reference,
                            const SE3& P_opt,
                            const Eigen::Matrix<double, 6, 6>& information,
                            const LiliOptions& options) {
    DetectorResult result;
    result.method = "lili";
    result.condition_number = 0.0;
    result.threshold = 0.0;

    double L = options.length_scale;
    if (L < 0.0) L = std::max(1e-6, medianRadius(local));

    // Algorithm 1, line 2: adaptively generated axis-aligned perturbations.
    const PerturbationScale scale =
        calibratePerturbation(local, reference, options.reassociation_k);
    const std::vector<Vector6d> dirs =
        generatePerturbations(scale, options.bidirectional);

    IcpOptions icp_opts = options.icp_options;
    if (icp_opts.max_iterations <= 0) icp_opts.max_iterations = 50;
    // The perturbation displaces points by roughly `scale.translation`; the
    // re-optimization can only "slide" if those points stay inside the
    // association gate. Enforce a floor tied to the perturbation so a stalled
    // (all-rejected) association cannot masquerade as a degenerate direction.
    const double assoc_floor = 2.0 * scale.translation;
    icp_opts.max_correspondence_distance =
        std::max(icp_opts.max_correspondence_distance, assoc_floor);

    // Algorithm 1, lines 3-7: re-optimize each perturbation and keep the twists
    // whose re-optimized displacement clears the tau_displacement gate.
    std::vector<Vector6d> D;
    D.reserve(dirs.size());
    for (const Vector6d& xi : dirs) {
        const double pert_disp = medianTwistDisplacement(xi, local);
        if (pert_disp < 1e-12) continue;

        const SE3 P_pert = P_opt * expSE3(xi);
        const IcpResult icp_res =
            pointToPlaneIcp(local, reference, P_pert, icp_opts);
        const SE3 T = P_opt.inverse() * icp_res.pose;
        const Vector6d t_deg = logSE3(T);
        const double reopt_disp = medianTwistDisplacement(t_deg, local);
        if (reopt_disp / pert_disp > options.tau_displacement) {
            D.push_back(t_deg);
        }
    }

    const int m = static_cast<int>(D.size());
    if (m == 0) {
        result.basis.resize(6, 0);
        result.eigenvalues.resize(0);
        result.eigenvectors.setZero();
        return result;
    }

    // Algorithm 1, line 8: PCA on the covariance of the detected directions
    // (Eq. 13). Rotation rows are weighted by the feature length first.
    Eigen::MatrixXd Rscaled(6, m);
    for (int i = 0; i < m; ++i) {
        Vector6d s = D[i];
        if (L != 0.0) s.head<3>() *= L;
        Rscaled.col(i) = s;
    }
    const Eigen::VectorXd mean = Rscaled.rowwise().mean();
    const Eigen::MatrixXd Rcent = Rscaled.colwise() - mean;
    const Eigen::MatrixXd Cov =
        (m > 1) ? Eigen::MatrixXd(Rcent * Rcent.transpose() /
                                  static_cast<double>(m - 1))
                : Eigen::MatrixXd::Zero(6, 6);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Cov);
    result.eigenvalues = solver.eigenvalues();  // ascending
    result.eigenvectors = solver.eigenvectors();

    // Algorithm 1, line 9: retain eigenvectors with eigenvalues ABOVE tau_PCA.
    const double lambda_max =
        result.eigenvalues.size() > 0 ? result.eigenvalues.maxCoeff() : 0.0;
    const double cut =
        std::max(options.pca_rel_threshold * lambda_max, options.pca_abs_threshold);
    result.threshold = cut;

    std::vector<int> deg_idx;
    for (int i = 0; i < result.eigenvalues.size(); ++i) {
        if (result.eigenvalues(i) > cut) deg_idx.push_back(i);
    }
    result.degenerate_indices = deg_idx;

    if (!deg_idx.empty()) {
        Eigen::Matrix<double, 6, Eigen::Dynamic> U(
            6, static_cast<int>(deg_idx.size()));
        for (size_t i = 0; i < deg_idx.size(); ++i) {
            U.col(static_cast<int>(i)) = result.eigenvectors.col(deg_idx[i]);
        }
        result.basis = l1SparsifyBasis(U, L, options.l1_iterations);
    } else {
        result.basis.resize(6, 0);
    }
    return result;
}

}  // namespace lili
