// M1 acceptance check: for every synthetic scene, the extended scan built from
// the analytically known degenerate basis must stay aligned (Q close to the
// plain-alignment baseline), while a control direction outside that basis must
// visibly degrade Q.

#include <Eigen/Dense>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "eval_quality.h"
#include "scenes.h"

using lili::PointCloud;
using lili::Scene;
using lili::SceneOptions;
using lili::Vector6d;

namespace {

// One direction orthogonal to the true degeneracy subspace, used as a control.
Eigen::Matrix<double, 6, Eigen::Dynamic> controlBasis(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& truth) {
    const Eigen::Matrix<double, 6, 6> I = Eigen::Matrix<double, 6, 6>::Identity();
    Eigen::Matrix<double, 6, Eigen::Dynamic> full;
    if (truth.cols() == 0) {
        full = I;
    } else {
        const Eigen::Matrix<double, 6, 6> Q =
            truth.householderQr().householderQ();
        full = Q.rightCols(6 - static_cast<int>(truth.cols()));
    }
    return full.col(0);
}

int reportScene(const Scene& scene, const SceneOptions& opt) {
    const auto q_none = lili::alignmentQualitySingle(scene.local, scene.reference, scene.P_gt);
    const auto q_true = lili::alignmentQuality(scene.local, scene.reference, scene.P_gt,
                                               scene.true_basis);
    const auto q_init = lili::alignmentQualitySingle(scene.local, scene.reference, scene.P_init);
    const auto control = controlBasis(scene.true_basis);
    const auto q_control = lili::alignmentQuality(scene.local, scene.reference, scene.P_gt, control);

    std::cout << "=== " << scene.name << " ===\n";
    std::cout << "  points: reference=" << scene.reference.size()
              << " local=" << scene.local.size() << "\n";
    std::cout << "  mean NN spacing (reference): "
              << lili::meanNearestNeighborSpacing(scene.reference) << "\n";
    std::cout << "  expected degeneracy: " << scene.expected_comment
              << " [dim=" << scene.true_basis.cols()
              << (scene.ground_truth_confirmed ? ", analytic" : ", hypothesis")
              << "]\n";
    std::cout << "  Q(P_init, unextended) = " << q_init.median << "\n";
    std::cout << "  Q(P_gt,   unextended) = " << q_none.median << "\n";
    std::cout << "  Q(true basis)         = " << q_true.median << "\n";
    std::cout << "  Q(control direction)  = " << q_control.median << "\n";
    // Two independent checks. The floor is the plain-alignment error or half a
    // point spacing, whichever is larger: a discrete reference cloud cannot be
    // hit exactly when the scan slides along the surface, so "stays on the
    // surface" means within sampling resolution, not exactly zero.
    const double spacing = lili::meanNearestNeighborSpacing(scene.reference);
    const double floor = std::max(q_none.median, 0.5 * spacing);
    const bool stays_on_surface = q_true.median <= 1.5 * floor;
    const bool separation = q_control.median > 1.5 * q_true.median;
    std::cout << "  stays_on_surface (Q_true <= 1.5*max(floor, h/2)): "
              << (stays_on_surface ? "PASS" : "FAIL") << "\n";
    std::cout << "  separation       (Q_control > 1.5*Q_true):        "
              << (separation ? "PASS" : "FAIL") << "\n";
    if (!scene.notes.empty()) { std::cout << "  note: " << scene.notes << "\n"; }
    std::cout << "\n";
    return (stays_on_surface && separation) ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    SceneOptions opt;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--noise" && i + 1 < argc) {
            opt.noise_sigma = std::atof(argv[++i]);
        } else if (arg == "--spacing" && i + 1 < argc) {
            opt.spacing = std::atof(argv[++i]);
        } else {
            std::cerr << "Usage: " << argv[0] << " [--noise sigma] [--spacing h]\n";
            return 1;
        }
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Scene report (noise sigma = " << opt.noise_sigma
              << ", spacing = " << opt.spacing << ")\n\n";

    int failures = 0;
    for (const Scene& scene : lili::makeAllScenes(opt)) {
        failures += reportScene(scene, opt);
    }
    std::cout << (failures == 0 ? "Scene acceptance: PASS"
                                : "Scene acceptance: FAIL")
              << " (" << failures << " failing scene(s))\n";
    return failures == 0 ? 0 : 1;
}
