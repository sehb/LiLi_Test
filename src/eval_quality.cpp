#include "eval_quality.h"

#include <algorithm>
#include <cmath>

namespace lili {
namespace {

constexpr double kPi = 3.14159265358979323846;

// van der Corput radical inverse in `base`.
double radicalInverse(int index, int base) {
    double f = 1.0, r = 0.0;
    int i = index + 1;
    while (i > 0) {
        f /= static_cast<double>(base);
        r += f * static_cast<double>(i % base);
        i /= base;
    }
    return r;
}

// Halton points in [0,1)^d, d <= 6, using the first six primes.
Eigen::VectorXd halton(int index, int dim) {
    static const int kPrimes[6] = {2, 3, 5, 7, 11, 13};
    Eigen::VectorXd c(dim);
    for (int j = 0; j < dim; ++j) {
        c(j) = 2.0 * radicalInverse(index, kPrimes[j % 6]) - 1.0;  // [-1, 1)
    }
    return c;
}

void collectNearestDistances(const PointCloud& query, const GridIndex& index,
                             double max_radius, std::vector<double>& out) {
    out.reserve(out.size() + query.size());
    for (const auto& p : query.points) {
        const double d = max_radius > 0.0 ? index.nearest(p, max_radius).distance
                                          : index.nearest(p).distance;
        // Censor at the search cap so Q stays bounded: beyond a few point
        // spacings every distance simply means "off the surface".
        out.push_back(std::isfinite(d) ? d : std::max(max_radius, 0.0));
    }
}

QualityResult summarize(std::vector<double>& d) {
    QualityResult r;
    if (d.empty()) { return r; }
    std::sort(d.begin(), d.end());
    double sum = 0.0;
    for (const double v : d) { sum += v; }
    r.median = d[d.size() / 2];
    r.mean = sum / static_cast<double>(d.size());
    r.p90 = d[std::min(d.size() - 1,
                       static_cast<std::size_t>(0.9 * static_cast<double>(d.size())))];
    r.samples = static_cast<int>(d.size());
    return r;
}

double boundingBoxScale(const PointCloud& reference) {
    if (reference.empty()) { return 1.0; }
    Eigen::Vector3d lo = reference.points.front(), hi = reference.points.front();
    for (const auto& p : reference.points) {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    return std::max(1e-6, (hi - lo).norm() / 64.0);
}

}  // namespace

QualityResult alignmentQuality(
    const PointCloud& local, const PointCloud& reference, const SE3& P_opt,
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& basis,
    const ExtensionOptions& opt) {
    if (local.empty() || reference.empty()) { return QualityResult{}; }

    const double spacing = medianNearestNeighborSpacing(reference, 200);
    const double radius =
        opt.max_radius > 0.0 ? opt.max_radius
                             : 8.0 * std::max(spacing, boundingBoxScale(reference) * 1e-3);
    // Sizing the cell to half the search radius bounds the ring scan to a few
    // rings even for points that end up far from the reference.
    const GridIndex index(reference, radius / 2.0);
    const int dim = static_cast<int>(basis.cols());
    std::vector<double> distances;
    distances.reserve(static_cast<std::size_t>(opt.n_samples + 1) * local.size());

    // Sample 0 is the unperturbed pose, so the extended scan always contains
    // the plain alignment as a baseline.
    for (int i = 0; i <= opt.n_samples; ++i) {
        Vector6d delta = Vector6d::Zero();
        if (dim > 0) {
            const Eigen::VectorXd c =
                opt.coefficient_scale * halton(i, dim);
            // Normalize basis columns to compare directions only (per handover.md)
            Eigen::Matrix<double, 6, Eigen::Dynamic> Bn = basis;
            for (int j = 0; j < Bn.cols(); ++j) {
                double n = Bn.col(j).norm();
                if (n > 1e-12) Bn.col(j) /= n;
            }
            delta = Bn * c;
        }
        const PointCloud extended = local.transformed(P_opt * expSE3(delta));
        collectNearestDistances(extended, index, radius, distances);
    }
    return summarize(distances);
}

QualityResult alignmentQualitySingle(const PointCloud& local,
                                     const PointCloud& reference,
                                     const SE3& T) {
    if (local.empty() || reference.empty()) { return QualityResult{}; }
    const double spacing = medianNearestNeighborSpacing(reference, 200);
    const double radius =
        8.0 * std::max(spacing, boundingBoxScale(reference) * 1e-3);
    const GridIndex index(reference, radius / 2.0);
    std::vector<double> distances;
    collectNearestDistances(local.transformed(T), index, radius, distances);
    return summarize(distances);
}

std::vector<double> principalAnglesDeg(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& A_in,
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& B_in, double length_scale) {
    std::vector<double> angles;
    if (A_in.cols() == 0 || B_in.cols() == 0) { return angles; }

    // Normalize columns in scaled space to compare directions only
    auto normalizeCols = [length_scale](Eigen::Matrix<double, 6, Eigen::Dynamic> M) {
        Eigen::Matrix<double, 6, Eigen::Dynamic> N = M;
        if (length_scale != 0.0) { N.topRows<3>() *= length_scale; }
        for (int j = 0; j < N.cols(); ++j) {
            double n = N.col(j).norm();
            if (n > 1e-12) N.col(j) /= n;
        }
        return N;
    };
    const Eigen::MatrixXd An = normalizeCols(A_in);
    const Eigen::MatrixXd Bn = normalizeCols(B_in);

    const Eigen::MatrixXd Aq = An.householderQr().householderQ();
    const Eigen::MatrixXd Bq = Bn.householderQr().householderQ();

    const int k = std::min<int>(static_cast<int>(A_in.cols()), static_cast<int>(B_in.cols()));
    if (k <= 0) return angles;
    const Eigen::MatrixXd M = Aq.leftCols(k).transpose() * Bq.leftCols(k);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M);
    for (int i = 0; i < k; ++i) {
        const double s = std::max(0.0, std::min(1.0, svd.singularValues()(i)));
        angles.push_back(std::acos(s) * 180.0 / kPi);
    }
    std::sort(angles.begin(), angles.end(), std::greater<double>());
    return angles;
}

SubspaceDeviation deviationFromSubspace(
    const Vector6d& twist,
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& basis_in,
    double length_scale) {
    SubspaceDeviation out;
    Vector6d scaled = twist;
    scaled.head<3>() *= length_scale;
    out.total = scaled.norm();
    if (out.total == 0.0) { return out; }

    if (basis_in.cols() == 0) {
        out.orthogonal = out.total;
        out.fraction = 1.0;
        return out;
    }

    Eigen::Matrix<double, 6, Eigen::Dynamic> scaled_basis = basis_in;
    if (length_scale != 0.0) { scaled_basis.topRows<3>() *= length_scale; }
    for (int j = 0; j < scaled_basis.cols(); ++j) {
        double n = scaled_basis.col(j).norm();
        if (n > 1e-12) scaled_basis.col(j) /= n;
    }
    Eigen::MatrixXd q_full = scaled_basis.householderQr().householderQ();
    Eigen::MatrixXd Q = q_full.leftCols(static_cast<int>(scaled_basis.cols()));
    const Vector6d projection = Q * (Q.transpose() * scaled);
    out.orthogonal = (scaled - projection).norm();
    out.fraction = out.orthogonal / out.total;
    return out;
}

}  // namespace lili
