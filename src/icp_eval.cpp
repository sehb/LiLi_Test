// M2 driver: register every scene with the shared ICP, then compare the
// Hessian-based (Zhang) detector against the analytically known degeneracy
// subspace using the paper's own metric Q.

#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "detectors.h"
#include "eval_quality.h"
#include "icp.h"
#include "scenes.h"

using lili::DetectorResult;
using lili::IcpOptions;
using lili::IcpResult;
using lili::Scene;
using lili::SceneOptions;
using lili::SE3;

namespace {

double rotationErrorDeg(const SE3& a, const SE3& b) {
    const SE3 e = a.inverse() * b;
    const double c = (e.R.trace() - 1.0) * 0.5;
    return std::acos(std::max(-1.0, std::min(1.0, c))) * 180.0 / M_PI;
}

double translationError(const SE3& a, const SE3& b) { return (a.t - b.t).norm(); }

double characteristicLength(const lili::PointCloud& local) {
    if (local.empty()) { return 1.0; }
    std::vector<double> d;
    d.reserve(local.size());
    for (const auto& p : local.points) { d.push_back(p.norm()); }
    std::sort(d.begin(), d.end());
    return std::max(1e-6, d[d.size() / 2]);
}

Eigen::Matrix<double, 6, Eigen::Dynamic> controlDirection(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& truth) {
    if (truth.cols() == 0) { return Eigen::Matrix<double, 6, 6>::Identity().col(0); }
    const Eigen::Matrix<double, 6, 6> Q = truth.householderQr().householderQ();
    return Q.rightCols(6 - static_cast<int>(truth.cols())).col(0);
}

struct SceneOutcome {
    std::string name;
    bool icp_ok = false;
    bool truth_dim_ok = false;
    int true_dim = 0;
    int zhang_dim = 0;
    double rot_err_deg = 0.0;
    double trans_err = 0.0;
    double err_total = 0.0;
    double err_orthogonal_fraction = 0.0;
    int iterations = 0;
    bool converged = false;
    double condition_number = 0.0;
    double q_true = 0.0;
    double q_zhang = 0.0;
    double q_control = 0.0;
    double q_plain = 0.0;
    double max_angle_deg = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
    SceneOptions scene_opt;
    IcpOptions icp_opt;
    int q_samples = 32;
    double pose_tol_trans = 0.01;
    double pose_tol_rot = 0.5;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--noise" && i + 1 < argc) {
            scene_opt.noise_sigma = std::atof(argv[++i]);
        } else if (arg == "--spacing" && i + 1 < argc) {
            scene_opt.spacing = std::atof(argv[++i]);
        } else if (arg == "--corr" && i + 1 < argc) {
            icp_opt.max_correspondence_distance = std::atof(argv[++i]);
        } else if (arg == "--q-samples" && i + 1 < argc) {
            q_samples = std::atoi(argv[++i]);
        } else if (arg == "--local-scale" && i + 1 < argc) {
            scene_opt.local_spacing_scale = std::atof(argv[++i]);
        } else if (arg == "--local-phase" && i + 1 < argc) {
            scene_opt.local_phase = std::atof(argv[++i]);
        } else {
            std::cerr << "Usage: " << argv[0]
                      << " [--noise s] [--spacing h] [--corr d] [--q-samples n]"
                      << " [--local-scale k] [--local-phase p]\n";
            return 1;
        }
    }

    lili::ExtensionOptions ext;
    ext.n_samples = q_samples;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "ICP + Zhang baseline (noise sigma = " << scene_opt.noise_sigma
              << ", spacing = " << scene_opt.spacing
              << ", max correspondence = " << icp_opt.max_correspondence_distance
              << ")\n";
    std::cout << "local view sampling: scale = " << scene_opt.local_spacing_scale
              << ", phase = " << scene_opt.local_phase
              << (scene_opt.local_spacing_scale == 1.0 && scene_opt.local_phase == 0.0
                      ? "  (exact subset of the reference)"
                      : "  (independently sampled)")
              << "\n\n";

    std::vector<SceneOutcome> outcomes;
    for (const Scene& scene : lili::makeAllScenes(scene_opt)) {
        SceneOutcome out;
        out.name = scene.name;
        out.true_dim = static_cast<int>(scene.true_basis.cols());

        const IcpResult icp =
            lili::pointToPlaneIcp(scene.local, scene.reference, scene.P_init, icp_opt);
        out.iterations = icp.iterations;
        out.converged = icp.converged;
        out.rot_err_deg = rotationErrorDeg(icp.pose, scene.P_gt);
        out.trans_err = translationError(icp.pose, scene.P_gt);

        // A degenerate scene does not pin down the pose, so the meaningful
        // acceptance test is that the pose error lies in the known subspace
        // (or, for the non-degenerate control, that it is essentially zero).
        const double length_scale = characteristicLength(scene.local);
        const lili::Vector6d err = logSE3(scene.P_gt.inverse() * icp.pose);
        const auto dev =
            lili::deviationFromSubspace(err, scene.true_basis, length_scale);
        out.err_total = dev.total;
        out.err_orthogonal_fraction = dev.fraction;
        if (scene_opt.noise_sigma > 0.0) {
            // With noise the pose error spreads over every weakly constrained
            // direction, so only boundedness is checkable.
            out.icp_ok = icp.converged && dev.total < 0.25;
        } else {
            out.icp_ok = icp.converged &&
                         (scene.true_basis.cols() == 0
                              ? dev.total < pose_tol_trans
                              : dev.fraction < 0.10);
        }

        const DetectorResult zhang = lili::zhangDetector(icp.information);
        out.zhang_dim = static_cast<int>(zhang.basis.cols());
        out.condition_number = zhang.condition_number;
        out.truth_dim_ok = (out.zhang_dim == out.true_dim);

        const auto q_plain =
            lili::alignmentQualitySingle(scene.local, scene.reference, icp.pose);
        const auto q_true = lili::alignmentQuality(scene.local, scene.reference,
                                                   icp.pose, scene.true_basis, ext);
        const auto q_zhang = lili::alignmentQuality(scene.local, scene.reference,
                                                    icp.pose, zhang.basis, ext);
        const auto q_control = lili::alignmentQuality(
            scene.local, scene.reference, icp.pose, controlDirection(scene.true_basis), ext);
        out.q_plain = q_plain.median;
        out.q_true = q_true.median;
        out.q_zhang = q_zhang.median;
        out.q_control = q_control.median;

        const auto angles =
            lili::principalAnglesDeg(zhang.basis, scene.true_basis, length_scale);
        out.max_angle_deg = angles.empty() ? 90.0 : angles.front();
        if (out.true_dim == 0) { out.max_angle_deg = 0.0; }

        std::cout << "=== " << scene.name << " ===\n";
        std::cout << "  ICP: iter=" << out.iterations
                  << " converged=" << (out.converged ? "yes" : "no")
                  << " corr=" << icp.correspondences << " rms=" << icp.rms << "\n";
        std::cout << "  pose error vs P_gt: rot=" << out.rot_err_deg
                  << " deg, trans=" << out.trans_err << " m\n";
        std::cout << "  pose error twist: |err|=" << out.err_total
                  << "  fraction outside true subspace=" << out.err_orthogonal_fraction
                  << "\n";
        std::cout << "  expected dim=" << out.true_dim
                  << "  Zhang dim=" << out.zhang_dim
                  << "  cond=" << out.condition_number << "\n";
        std::cout << "  Q: plain=" << out.q_plain << " true=" << out.q_true
                  << " zhang=" << out.q_zhang << " control=" << out.q_control << "\n";
        std::cout << "  max principal angle (Zhang vs truth) = " << out.max_angle_deg
                  << " deg\n";
        std::cout << "  ICP " << (out.icp_ok ? "OK" : "FAILED")
                  << " | dim match " << (out.truth_dim_ok ? "yes" : "no") << "\n\n";

        outcomes.push_back(out);
    }

    int icp_failures = 0, dim_mismatch = 0;
    for (const auto& o : outcomes) {
        if (!o.icp_ok) { ++icp_failures; }
        if (!o.truth_dim_ok) { ++dim_mismatch; }
    }
    std::cout << "M2 acceptance: ICP " << (icp_failures == 0 ? "PASS" : "FAIL")
              << " (" << icp_failures << " scene(s) failed to register)\n";
    std::cout << "Zhang dimension match: " << (outcomes.size() - dim_mismatch) << "/"
              << outcomes.size() << " scenes (mismatches are the finding, not a test failure)\n";
    return icp_failures == 0 ? 0 : 1;
}
