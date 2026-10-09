#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <vector>

#include "se3.h"

namespace lili {

struct PointCloud {
    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector3d> normals;  // optional; parallel to `points` when set

    std::size_t size() const { return points.size(); }
    bool empty() const { return points.empty(); }
    bool hasNormals() const { return normals.size() == points.size(); }

    void transform(const SE3& T);
    PointCloud transformed(const SE3& T) const;
    void append(const PointCloud& other);
    void reserve(std::size_t n);
};

struct NeighborResult {
    int index = -1;
    double distance = std::numeric_limits<double>::infinity();
};

// Uniform-grid nearest-neighbour index. O(1)-ish queries, adequate for the
// point counts used by the synthetic scenes (no PCL dependency by design).
class GridIndex {
public:
    explicit GridIndex(const PointCloud& cloud, double cell_size);

    // Returns the closest point in `cloud`, expanding rings until the
    // remaining cells provably cannot hold a closer point.
    NeighborResult nearest(const Eigen::Vector3d& query,
                           double max_radius =
                               std::numeric_limits<double>::infinity()) const;

    double cellSize() const { return cell_size_; }

private:
    struct CellKey {
        int x = 0, y = 0, z = 0;
        bool operator==(const CellKey& o) const {
            return x == o.x && y == o.y && z == o.z;
        }
    };
    struct CellHash {
        std::size_t operator()(const CellKey& k) const;
    };

    CellKey cellOf(const Eigen::Vector3d& p) const;

    const PointCloud& cloud_;
    double cell_size_;
    Eigen::Vector3d origin_ = Eigen::Vector3d::Zero();
    std::unordered_map<CellKey, std::vector<int>, CellHash> cells_;
};

// Median/mean nearest-neighbour spacing, the diagnostic the paper uses to
// characterise its synthetic clouds (0.05 in the reference setup).
double meanNearestNeighborSpacing(const PointCloud& cloud, int max_samples = 300);
double medianNearestNeighborSpacing(const PointCloud& cloud, int max_samples = 300);

}  // namespace lili
