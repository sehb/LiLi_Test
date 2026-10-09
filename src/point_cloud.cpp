#include "point_cloud.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace lili {

void PointCloud::transform(const SE3& T) {
    for (auto& p : points) { p = T * p; }
    for (auto& n : normals) { n = T.R * n; }
}

PointCloud PointCloud::transformed(const SE3& T) const {
    PointCloud out = *this;
    out.transform(T);
    return out;
}

void PointCloud::append(const PointCloud& other) {
    points.insert(points.end(), other.points.begin(), other.points.end());
    const bool both_have_normals = hasNormals() && other.hasNormals();
    if (both_have_normals) {
        normals.insert(normals.end(), other.normals.begin(), other.normals.end());
    } else {
        normals.clear();
    }
}

void PointCloud::reserve(std::size_t n) {
    points.reserve(n);
    normals.reserve(n);
}

std::size_t GridIndex::CellHash::operator()(const CellKey& k) const {
    // 64-bit mix of the three cell coordinates.
    std::uint64_t h = 1469598103934665603ULL;
    const auto mix = [&h](int v) {
        h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(v));
        h *= 1099511628211ULL;
    };
    mix(k.x);
    mix(k.y);
    mix(k.z);
    return static_cast<std::size_t>(h);
}

GridIndex::GridIndex(const PointCloud& cloud, double cell_size)
    : cloud_(cloud), cell_size_(cell_size > 0.0 ? cell_size : 1.0) {
    if (cloud_.points.empty()) { return; }

    Eigen::Vector3d lo = cloud_.points.front();
    for (const auto& p : cloud_.points) { lo = lo.cwiseMin(p); }
    origin_ = lo;

    for (std::size_t i = 0; i < cloud_.points.size(); ++i) {
        cells_[cellOf(cloud_.points[i])].push_back(static_cast<int>(i));
    }
}

GridIndex::CellKey GridIndex::cellOf(const Eigen::Vector3d& p) const {
    CellKey k;
    k.x = static_cast<int>(std::floor((p.x() - origin_.x()) / cell_size_));
    k.y = static_cast<int>(std::floor((p.y() - origin_.y()) / cell_size_));
    k.z = static_cast<int>(std::floor((p.z() - origin_.z()) / cell_size_));
    return k;
}

NeighborResult GridIndex::nearest(const Eigen::Vector3d& query,
                                  double max_radius) const {
    NeighborResult result;
    if (cloud_.points.empty()) { return result; }

    const double max_r2 =
        std::isfinite(max_radius) ? max_radius * max_radius : -1.0;
    const CellKey center = cellOf(query);

    // Ring 0 is the query cell itself; expand outward and stop once the
    // closest remaining cell boundary is farther than the best hit.
    const int max_ring = std::isfinite(max_radius)
                             ? static_cast<int>(std::ceil(max_radius / cell_size_)) + 1
                             : 1000;

    for (int r = 0; r <= max_ring; ++r) {
        for (int dx = -r; dx <= r; ++dx) {
            for (int dy = -r; dy <= r; ++dy) {
                for (int dz = -r; dz <= r; ++dz) {
                    if (std::max(std::abs(dx), std::max(std::abs(dy), std::abs(dz))) != r) {
                        continue;
                    }
                    const CellKey key{center.x + dx, center.y + dy, center.z + dz};
                    const auto it = cells_.find(key);
                    if (it == cells_.end()) { continue; }
                    for (const int idx : it->second) {
                        const double d2 =
                            (cloud_.points[idx] - query).squaredNorm();
                        if (max_r2 >= 0.0 && d2 > max_r2) { continue; }
                        if (d2 < result.distance * result.distance) {
                            result.distance = std::sqrt(d2);
                            result.index = idx;
                        }
                    }
                }
            }
        }
        // Points in ring r+1 are at least r * cell_size away.
        if (result.index >= 0 && (static_cast<double>(r) * cell_size_) > result.distance) {
            break;
        }
    }
    return result;
}

namespace {

// Deterministic thinning so diagnostics do not depend on a random seed.
int strideFor(std::size_t n, int max_samples) {
    if (max_samples <= 0 || n <= static_cast<std::size_t>(max_samples)) { return 1; }
    return static_cast<int>(n / static_cast<std::size_t>(max_samples));
}

double spacingStatistic(const PointCloud& cloud, int max_samples, bool median) {
    if (cloud.size() < 2) { return 0.0; }
    const int stride = strideFor(cloud.size(), max_samples);

    std::vector<double> d;
    d.reserve(cloud.size() / static_cast<std::size_t>(stride) + 1);
    for (std::size_t i = 0; i < cloud.size();
         i += static_cast<std::size_t>(stride)) {
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < cloud.size(); ++j) {
            if (j == i) { continue; }
            const double dist = (cloud.points[j] - cloud.points[i]).squaredNorm();
            if (dist > 0.0 && dist < best) { best = dist; }
        }
        if (std::isfinite(best)) { d.push_back(std::sqrt(best)); }
    }
    if (d.empty()) { return 0.0; }
    if (median) {
        std::sort(d.begin(), d.end());
        return d[d.size() / 2];
    }
    double sum = 0.0;
    for (const double v : d) { sum += v; }
    return sum / static_cast<double>(d.size());
}

}  // namespace

double meanNearestNeighborSpacing(const PointCloud& cloud, int max_samples) {
    return spacingStatistic(cloud, max_samples, false);
}

double medianNearestNeighborSpacing(const PointCloud& cloud, int max_samples) {
    return spacingStatistic(cloud, max_samples, true);
}

}  // namespace lili
