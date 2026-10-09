// M3 driver: register every scene with the shared ICP, then compare three
// detectors against the analytically known degeneracy subspace using the
// paper's own metric Q:
//   - the Hessian-based Zhang baseline,
//   - the perturbation-based LiLi detector,
//   - LiLi with frozen correspondences (the H3 re-association ablation).

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "detectors.h"
#include "eval_quality.h"
#include "icp.h"
#include "lili_detector.h"
#include "scenes.h"

using lili::DetectorResult;
using lili::IcpOptions;
using lili::IcpResult;
using lili::LiliOptions;
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

// Max principal angle in degrees between a detected basis and the truth, using
// the scene's characteristic length; 0 when the scene is non-degenerate.
double maxAngleDeg(const Eigen::Matrix<double, 6, Eigen::Dynamic>& basis,
                   const Eigen::Matrix<double, 6, Eigen::Dynamic>& truth,
                   double length_scale) {
    if (truth.cols() == 0 || basis.cols() == 0) { return 0.0; }
    const std::vector<double> angles =
        lili::principalAnglesDeg(basis, truth, length_scale);
    return angles.empty() ? 90.0 : angles.front();
}

struct Row {
    std::string name;
    int true_dim = 0;
    bool icp_ok = false;
    int zhang_dim = 0;
    int lili_dim = 0;
    int ablated_dim = 0;
    double zhang_angle = 0.0;
    double lili_angle = 0.0;
    double ablated_angle = 0.0;
    double q_true = 0.0;
    double q_zhang = 0.0;
    double q_lili = 0.0;
    double q_ablated = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
    SceneOptions scene_opt;
    IcpOptions icp_opt;
    LiliOptions lili_opt;
    int q_samples = 32;

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
        } else if (arg == "--k" && i + 1 < argc) {
            lili_opt.reassociation_k = std::atoi(argv[++i]);
        } else if (arg == "--tau-disp" && i + 1 < argc) {
            lili_opt.tau_displacement = std::atof(argv[++i]);
        } else if (arg == "--pca-rel" && i + 1 < argc) {
            lili_opt.pca_rel_threshold = std::atof(argv[++i]);
        } else {
            std::cerr << "Usage: " << argv[0]
                      << " [--noise s] [--spacing h] [--corr d] [--q-samples n]"
                      << " [--local-scale k] [--local-phase p] [--k k]"
                      << " [--tau-disp t] [--pca-rel r]\n";
            return 1;
        }
    }

    lili::ExtensionOptions ext;
    ext.n_samples = q_samples;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "M3 comparison (noise sigma = " << scene_opt.noise_sigma
              << ", spacing = " << scene_opt.spacing << ", k = "
              << lili_opt.reassociation_k << ")\n";
    std::cout << "local view sampling: scale = " << scene_opt.local_spacing_scale
              << ", phase = " << scene_opt.local_phase
              << (scene_opt.local_spacing_scale == 1.0 && scene_opt.local_phase == 0.0
                      ? "  (exact subset of the reference)"
                      : "  (independently sampled)")
              << "\n\n";

    std::vector<Row> rows;
    int icp_failures = 0;
    for (const Scene& scene : lili::makeAllScenes(scene_opt)) {
        Row row;
        row.name = scene.name;
        row.true_dim = static_cast<int>(scene.true_basis.cols());

        const IcpResult icp =
            lili::pointToPlaneIcp(scene.local, scene.reference, scene.P_init, icp_opt);
        const double length_scale = characteristicLength(scene.local);
        const lili::Vector6d err = logSE3(scene.P_gt.inverse() * icp.pose);
        const auto dev =
            lili::deviationFromSubspace(err, scene.true_basis, length_scale);
        if (scene_opt.noise_sigma > 0.0) {
            row.icp_ok = icp.converged && dev.total < 0.25;
        } else {
            row.icp_ok = icp.converged &&
                         (scene.true_basis.cols() == 0
                              ? dev.total < 0.01
                              : dev.fraction < 0.10);
        }
        if (!row.icp_ok) { ++icp_failures; }

        // Detectors.
        const DetectorResult zhang = lili::zhangDetector(icp.information);
        const DetectorResult li = lili::liliDetector(
            scene.local, scene.reference, icp.pose, icp.information, lili_opt);
        LiliOptions ablated = lili_opt;
        ablated.icp_options.reassociate = false;
        DetectorResult li_abl = lili::liliDetector(
            scene.local, scene.reference, icp.pose, icp.information, ablated);

        row.zhang_dim = static_cast<int>(zhang.basis.cols());
        row.lili_dim = static_cast<int>(li.basis.cols());
        row.ablated_dim = static_cast<int>(li_abl.basis.cols());
        row.zhang_angle = maxAngleDeg(zhang.basis, scene.true_basis, length_scale);
        row.lili_angle = maxAngleDeg(li.basis, scene.true_basis, length_scale);
        row.ablated_angle = maxAngleDeg(li_abl.basis, scene.true_basis, length_scale);

        row.q_true = lili::alignmentQuality(scene.local, scene.reference, icp.pose,
                                            scene.true_basis, ext)
                         .median;
        row.q_zhang = lili::alignmentQuality(scene.local, scene.reference, icp.pose,
                                             zhang.basis, ext)
                          .median;
        row.q_lili = lili::alignmentQuality(scene.local, scene.reference, icp.pose,
                                            li.basis, ext)
                         .median;
        row.q_ablated = lili::alignmentQuality(scene.local, scene.reference, icp.pose,
                                               li_abl.basis, ext)
                            .median;

        std::cout << "=== " << scene.name << " ===\n";
        std::cout << "  dims: true=" << row.true_dim << " zhang=" << row.zhang_dim
                  << " lili=" << row.lili_dim << " lili(ablated)=" << row.ablated_dim
                  << "\n";
        std::cout << "  max principal angle: zhang=" << row.zhang_angle
                  << " lili=" << row.lili_angle << " ablated=" << row.ablated_angle
                  << " deg\n";
        std::cout << "  Q: true=" << row.q_true << " zhang=" << row.q_zhang
                  << " lili=" << row.q_lili << " ablated=" << row.q_ablated << "\n\n";
        rows.push_back(row);
    }

    // Compact comparison table + H2-style summary.
    std::cout << std::setw(20) << "scene" << std::setw(6) << "tru" << std::setw(6)
              << "zD" << std::setw(6) << "lD" << std::setw(8) << "zAng" << std::setw(8)
              << "lAng" << std::setw(9) << "Qzhang" << std::setw(9) << "Qlili"
              << std::setw(10) << "Qablated" << "\n";
    int lili_better = 0;
    for (const Row& r : rows) {
        std::cout << std::setw(20) << r.name << std::setw(6) << r.true_dim
                  << std::setw(6) << r.zhang_dim << std::setw(6) << r.lili_dim
                  << std::setw(8) << r.zhang_angle << std::setw(8) << r.lili_angle
                  << std::setw(9) << r.q_zhang << std::setw(9) << r.q_lili
                  << std::setw(10) << r.q_ablated << "\n";
        if (r.true_dim > 0 && r.q_lili <= 0.5 * r.q_zhang) { ++lili_better; }
    }

    std::cout << "\nICP " << (icp_failures == 0 ? "PASS" : "FAIL") << " ("
              << icp_failures << " scene(s) failed to register)\n";
    if (scene_opt.noise_sigma == 0.0) {
        std::cout << "H1 (noise-free angle <= 5 deg, LiLi): "
                  << [&] {
                       for (const Row& r : rows) {
                         if (r.true_dim > 0 && r.lili_angle > 5.0) { return "no"; }
                       }
                       return "yes";
                     }()
                  << "\n";
    }
    std::cout << "H2 (Q_lili <= 0.5 Q_zhang): " << lili_better << "/" << rows.size()
              << " scenes\n";
    return icp_failures == 0 ? 0 : 1;
}
