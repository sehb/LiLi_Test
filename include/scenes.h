#pragma once

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "point_cloud.h"
#include "se3.h"

namespace lili {

// A synthetic registration problem with a known degeneracy subspace.
//
// Conventions:
//  - Shapes are generated in a "shape frame" whose z axis is the shape's
//    symmetry axis, then rotated by a fixed misaligning rotation so that no
//    scene axis coincides with x/y/z (as the paper deliberately does).
//  - `local` is the partial view expressed in the sensor frame; `reference`
//    is the full structure in the world frame.
//  - `P_gt` maps local -> reference exactly.
//  - `true_basis` columns live in se(3) expressed in the *body* (sensor)
//    frame, matching the LiLi basis convention of Eq. (7)/(9).
struct Scene {
    std::string name;
    PointCloud reference;
    PointCloud local;
    SE3 P_gt;
    SE3 P_init;                                     // starting guess for the ICP
    Eigen::Matrix<double, 6, Eigen::Dynamic> true_basis;
    bool ground_truth_confirmed = true;             // false = hypothesis only
    std::string notes;
    std::string expected_comment;                   // e.g. "2 translations + 1 rotation"
};

struct SceneOptions {
    double spacing = 0.05;   // target surface sampling step (paper uses 0.05)
    double noise_sigma = 0.0;
    int max_points = 8000;   // safety cap on the shared surface sampling
    unsigned seed = 20261009u;
    bool apply_initial_offset = true;  // P_init = exp(delta0) * P_gt
    // Local-view sampling. The defaults (1.0, 0.0) make the local view an exact
    // subset of the reference, which is the idealised case. Real scans are
    // independently sampled, so setting a different spacing and/or a half-step
    // phase makes the nearest-neighbour associations systematically offset -
    // the data-reassociation ambiguity the LiLi paper is about.
    double local_spacing_scale = 1.0;
    double local_phase = 0.0;  // in steps of the local spacing
};

Scene makePlaneScene(const SceneOptions& opt);
Scene makeClosedCylinderScene(const SceneOptions& opt);
Scene makeOpenCylinderScene(const SceneOptions& opt);
Scene makeSinusoidalCylinderScene(const SceneOptions& opt);
// Coupled-degeneracy scenario: a screw-symmetric tube, standing in for the
// paper's Rotating Tunnel. The exact symmetry is a single screw twist.
Scene makeScrewTubeScene(const SceneOptions& opt);
// Non-degenerate control: a rough height field with no symmetry (replaces the
// paper's Random3D for the false-positive check).
Scene makeRandomScene(const SceneOptions& opt);

std::vector<Scene> makeAllScenes(const SceneOptions& opt);

// Fixed misaligning rotation shared by all scenes.
Eigen::Matrix3d shapeToWorldRotation();

}  // namespace lili
