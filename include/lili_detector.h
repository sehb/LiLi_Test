#pragma once

#include <Eigen/Dense>
#include <vector>

#include "detectors.h"
#include "icp.h"
#include "point_cloud.h"
#include "se3.h"

namespace lili {

// LiLi detector options, following paper §V / Algorithm 1.
struct LiliOptions {
    // Target re-association distance k (paper §V). A single integer that sets
    // how far the perturbation should push a benchmark point (one at the median
    // radius from the scan center) so that its nearest neighbour in the
    // reference slides by k samples. Drives the adaptive perturbation scale.
    int reassociation_k = 3;
    // PCA eigenvalue threshold (paper line 9): keep eigenvectors whose
    // eigenvalue is *greater* than tau_PCA. Expressed relative to the largest
    // eigenvalue so the (unit-mixed) se(3) twists need no absolute scale.
    double pca_rel_threshold = 3e-2;
    // Absolute floor applied on top of the relative threshold.
    double pca_abs_threshold = 0.0;
    // Displacement gate (paper lines 6-7): a re-optimised point displacement,
    // normalised by the perturbation displacement, above this ratio marks the
    // perturbation as degenerate. Perturbations that get pulled back by the
    // re-optimisation are discarded before PCA.
    double tau_displacement = 0.10;
    // Feature length (m per rad) used to weigh the rotational rows of a twist
    // against the translational rows inside PCA and sparsification. 1.0 keeps
    // the raw se(3) coordinates the paper works in; <0 selects the local
    // cloud's median radius.
    double length_scale = 1.0;
    // Jacobi sweeps for the Eq. (14) L1 sparsification (0 disables it).
    int l1_iterations = 50;
    // Apply each axis couple with both signs (6 perturbations). The paper's set
    // is 3 couples; the sign-pair is what lets a 3-dimensional plane degeneracy
    // (needing >= 4 samples after mean removal) be observable.
    bool bidirectional = true;
    // ICP options used for the per-perturbation re-optimization.
    IcpOptions icp_options;
};

// LiLi detector: axis-aligned coupled perturbations -> per-perturbation
// re-optimization -> T = P_opt^-1 P_perturbed -> log -> tau_displacement gate ->
// PCA (Eq. 13) -> keep eigenvalues above tau_PCA -> L1 sparsification (Eq. 14).
DetectorResult liliDetector(const PointCloud& local,
                            const PointCloud& reference,
                            const SE3& P_opt,
                            const Eigen::Matrix<double, 6, 6>& information,
                            const LiliOptions& options = LiliOptions{});

// Adaptive perturbation magnitude for one axis couple (paper §V): both a
// translation (m) and a rotation (rad) large enough to slide a benchmark point
// at the local cloud's median radius by k reference-scan samples.
struct PerturbationScale {
    double translation = 0.0;  // metres along the axis
    double rotation = 0.0;     // radians about the axis
    double spacing = 0.0;      // median NN spacing of the reference scan
    double benchmark_radius = 0.0;
};

PerturbationScale calibratePerturbation(const PointCloud& local,
                                        const PointCloud& reference,
                                        int k);

// Axis-aligned perturbations (∆t_x,∆r_x),(∆t_y,∆r_y),(∆t_z,∆r_z) in body-frame
// se(3) (xi = [phi; rho]), scaled by `scale`. When `bidirectional` is set each
// couple is emitted with both signs.
std::vector<Vector6d> generatePerturbations(const PerturbationScale& scale,
                                            bool bidirectional);

// Eq. (14): minimise ||T B||_1 subject to T^T T = I, returning B_sparse = T B.
// Implemented with Givens (Jacobi) sweeps, so every step is an orthogonal
// transformation of the basis and the spanned subspace is preserved exactly.
Eigen::Matrix<double, 6, Eigen::Dynamic> l1SparsifyBasis(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& B,
    double length_scale,
    int max_iter);

}  // namespace lili
