#include "icp.h"

#include <algorithm>
#include <cmath>

namespace lili {
namespace {

double gridCellFor(const PointCloud& cloud, double max_radius) {
    if (cloud.empty()) { return 1.0; }
    Eigen::Vector3d lo = cloud.points.front(), hi = cloud.points.front();
    for (const auto& p : cloud.points) {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    const double cell = (max_radius > 0.0) ? 0.5 * max_radius
                                           : (hi - lo).norm() / 128.0;
    return std::max(cell, 1e-6);
}

struct Association {
    std::vector<int> local_indices;
    std::vector<int> reference_indices;
    int rejected_by_distance = 0;
    int rejected_by_normal = 0;
};

Association buildAssociation(const PointCloud& local, const PointCloud& reference,
                             const GridIndex& index, const SE3& pose,
                             const IcpOptions& options) {
    Association assoc;
    const bool use_normals =
        reference.hasNormals() && local.hasNormals() && options.normal_agreement > 0.0;
    const double max_dist = options.max_correspondence_distance;

    for (std::size_t i = 0; i < local.size(); ++i) {
        const Eigen::Vector3d p_world = pose * local.points[i];
        const NeighborResult nb = index.nearest(p_world, max_dist);
        if (nb.index < 0) {
            ++assoc.rejected_by_distance;
            continue;
        }
        if (use_normals) {
            const double agreement =
                reference.normals[nb.index].dot(pose.R * local.normals[i]);
            if (agreement < options.normal_agreement) {
                ++assoc.rejected_by_normal;
                continue;
            }
        }
        assoc.local_indices.push_back(static_cast<int>(i));
        assoc.reference_indices.push_back(nb.index);
    }
    return assoc;
}

// Accumulate the body-frame Gauss-Newton system at a fixed pose.
void accumulate(const PointCloud& local, const PointCloud& reference,
                const Association& assoc, const SE3& pose,
                Eigen::Matrix<double, 6, 6>* H, Eigen::Matrix<double, 6, 1>* g,
                double* sse) {
    H->setZero();
    g->setZero();
    *sse = 0.0;
    for (std::size_t k = 0; k < assoc.local_indices.size(); ++k) {
        const int i = assoc.local_indices[k];
        const int j = assoc.reference_indices[k];
        const Eigen::Vector3d& p = local.points[i];
        const Eigen::Vector3d& q = reference.points[j];
        const Eigen::Vector3d& n = reference.normals[j];
        const Eigen::Vector3d p_world = pose * p;

        const double r = n.dot(p_world - q);
        Eigen::Matrix<double, 1, 6> J;
        J.head<3>() = (-n.transpose() * pose.R * skew(p));
        J.tail<3>() = n.transpose() * pose.R;
        *H += J.transpose() * J;
        *g += J.transpose() * r;
        *sse += r * r;
    }
}

// One association + accumulation pass; returns the mean squared residual, which
// is the quantity the Levenberg-Marquardt accept/reject test compares.
bool evaluateAt(const PointCloud& local, const PointCloud& reference,
                const GridIndex& index, const SE3& pose, const IcpOptions& options,
                Eigen::Matrix<double, 6, 6>* H, Eigen::Matrix<double, 6, 1>* g,
                double* cost, int* count) {
    const Association assoc = buildAssociation(local, reference, index, pose, options);
    *count = static_cast<int>(assoc.local_indices.size());
    if (*count < 6) {
        H->setZero();
        g->setZero();
        *cost = std::numeric_limits<double>::infinity();
        return false;
    }
    double sse = 0.0;
    accumulate(local, reference, assoc, pose, H, g, &sse);
    *cost = sse / static_cast<double>(*count);
    return true;
}

double meanDiagonal(const Eigen::Matrix<double, 6, 6>& H) {
    return H.diagonal().mean();
}

}  // namespace

IcpResult pointToPlaneIcp(const PointCloud& local, const PointCloud& reference,
                          const SE3& initial, const IcpOptions& options) {
    IcpResult result;
    result.pose = initial;
    if (local.empty() || reference.empty() || !reference.hasNormals()) {
        return result;
    }

    const GridIndex index(reference,
                          gridCellFor(reference, options.max_correspondence_distance));

    result.pose = initial;
    double lambda = -1.0;

    for (int iter = 0; iter < options.max_iterations; ++iter) {
        Eigen::Matrix<double, 6, 6> H;
        Eigen::Matrix<double, 6, 1> g;
        double cost = 0.0;
        int count = 0;
        if (!evaluateAt(local, reference, index, result.pose, options, &H, &g, &cost,
                        &count)) {
            result.iterations = iter;
            break;  // too few correspondences to constrain anything
        }
        if (lambda < 0.0) {
            // Scale the initial damping to the problem so it is unit-agnostic.
            lambda = 1e-3 * std::max(1e-12, meanDiagonal(H));
        }
        result.rms = std::sqrt(cost);
        result.rms_history.push_back(result.rms);
        result.iterations = iter + 1;

        // Undamped Gauss-Newton on a degenerate scene produces unbounded steps
        // along the null space, so every step is damped and only accepted when
        // it actually reduces the mean squared residual.
        bool accepted = false;
        Eigen::Matrix<double, 6, 1> delta = Eigen::Matrix<double, 6, 1>::Zero();
        for (int trial = 0; trial < 12; ++trial) {
            const Eigen::Matrix<double, 6, 6> damped =
                H + lambda * Eigen::Matrix<double, 6, 6>::Identity();
            delta = damped.ldlt().solve(-g);
            if (!delta.allFinite()) { lambda *= 10.0; continue; }

            const SE3 candidate = result.pose * expSE3(delta);
            Eigen::Matrix<double, 6, 6> Hc;
            Eigen::Matrix<double, 6, 1> gc;
            double cost_c = 0.0;
            int count_c = 0;
            if (!evaluateAt(local, reference, index, candidate, options, &Hc, &gc,
                            &cost_c, &count_c) ||
                cost_c > cost) {
                lambda *= 10.0;
                continue;
            }
            result.pose = candidate;
            lambda = std::max(1e-15, lambda * 0.3);
            accepted = true;
            break;
        }

        if (!accepted) {
            result.converged = true;  // no damped step improves the cost any more
            break;
        }
        if (delta.head<3>().norm() < options.rotation_tolerance &&
            delta.tail<3>().norm() < options.translation_tolerance) {
            result.converged = true;
            break;
        }
    }

    // Re-evaluate at the final pose so the returned information matrix and
    // statistics describe exactly the returned transform.
    Eigen::Matrix<double, 6, 1> g;
    double cost = 0.0;
    int count = 0;
    if (evaluateAt(local, reference, index, result.pose, options, &result.information,
                   &g, &cost, &count)) {
        result.correspondences = count;
        result.rms = std::sqrt(cost);
    } else {
        result.correspondences = count;
        result.information.setZero();
    }
    return result;
}

AssociationStats associateAndEvaluate(const PointCloud& local,
                                      const PointCloud& reference,
                                      const SE3& pose,
                                      const IcpOptions& options) {
    AssociationStats stats;
    if (local.empty() || reference.empty() || !reference.hasNormals()) {
        return stats;
    }
    const GridIndex index(reference,
                          gridCellFor(reference, options.max_correspondence_distance));
    const Association assoc = buildAssociation(local, reference, index, pose, options);
    stats.correspondences = static_cast<int>(assoc.local_indices.size());
    stats.rejected_by_distance = assoc.rejected_by_distance;
    stats.rejected_by_normal = assoc.rejected_by_normal;

    double sse = 0.0;
    for (std::size_t k = 0; k < assoc.local_indices.size(); ++k) {
        const Eigen::Vector3d p_world = pose * local.points[assoc.local_indices[k]];
        const Eigen::Vector3d& q = reference.points[assoc.reference_indices[k]];
        const Eigen::Vector3d& n = reference.normals[assoc.reference_indices[k]];
        const double r = n.dot(p_world - q);
        sse += r * r;
    }
    if (stats.correspondences > 0) {
        stats.rms = std::sqrt(sse / static_cast<double>(stats.correspondences));
    }
    return stats;
}

}  // namespace lili
