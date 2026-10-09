#pragma once

#include <Eigen/Dense>
#include <vector>

#include "point_cloud.h"
#include "se3.h"

namespace lili {

// Point-to-plane Gauss-Newton ICP with re-association at every iteration.
//
// Convention: the increment acts on the right in the body frame,
//     T <- T * exp(delta^),   delta = [phi; rho]
// so the residual Jacobian of r = n^T (R p + t - q) is
//     J = [ -n^T R [p]_x , n^T R ].
// This makes the accumulated information matrix live in the *body* frame,
// which is exactly the frame used by the LiLi basis (paper Eq. 7 and 9).
//
// Note: this differs from the legacy InformationMatrixCalculator::addPointPlane,
// whose translation block is `n^T` (a world-frame translation increment mixed
// with a body-frame rotation). Do not treat the two as interchangeable.
struct IcpOptions {
    int max_iterations = 60;
    double max_correspondence_distance = 0.15;
    // Reject pairs whose normals disagree by more than this cosine (<=0 disables).
    double normal_agreement = 0.5;
    double rotation_tolerance = 1e-8;     // rad
    double translation_tolerance = 1e-8;  // m
};

struct IcpResult {
    SE3 pose;
    // Body-frame information matrix H = sum J^T J at the final correspondences.
    Eigen::Matrix<double, 6, 6> information = Eigen::Matrix<double, 6, 6>::Zero();
    int iterations = 0;
    int correspondences = 0;
    bool converged = false;
    double rms = 0.0;
    std::vector<double> rms_history;
};

// Requires `reference.normals` to be populated.
IcpResult pointToPlaneIcp(const PointCloud& local, const PointCloud& reference,
                          const SE3& initial,
                          const IcpOptions& options = IcpOptions{});

// Residual statistics at a fixed pose, using the same association rules.
struct AssociationStats {
    int correspondences = 0;
    int rejected_by_distance = 0;
    int rejected_by_normal = 0;
    double rms = 0.0;
};

AssociationStats associateAndEvaluate(const PointCloud& local,
                                      const PointCloud& reference,
                                      const SE3& pose,
                                      const IcpOptions& options = IcpOptions{});

}  // namespace lili
