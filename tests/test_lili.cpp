#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "eval_quality.h"
#include "icp.h"
#include "lili_detector.h"
#include "scenes.h"

namespace {

int failures = 0;

void expect(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    }
}

double maxAngle(const std::vector<double>& angles) {
    double m = 0.0;
    for (double a : angles) m = std::max(m, a);
    return m;
}

}  // namespace

int main() {
    using namespace lili;

    // 1. H1: on noise-free data the LiLi subspace must match the analytic truth
    //    (max principal angle <= 5 deg); the non-degenerate control must report
    //    no degeneracy (H4).
    {
        SceneOptions opt;
        opt.noise_sigma = 0.0;
        for (const Scene& scene : makeAllScenes(opt)) {
            const IcpResult icp =
                pointToPlaneIcp(scene.local, scene.reference, scene.P_init);
            const DetectorResult li = liliDetector(
                scene.local, scene.reference, icp.pose, icp.information);
            const int truth_dim = static_cast<int>(scene.true_basis.cols());
            const int lili_dim = static_cast<int>(li.basis.cols());
            std::cout << "  " << scene.name << ": truth_dim=" << truth_dim
                      << " lili_dim=" << lili_dim;
            if (truth_dim == 0) {
                std::cout << " (expect 0)\n";
                expect(lili_dim == 0, scene.name + ": no false positive (H4)");
                continue;
            }
            const double angle =
                maxAngle(principalAnglesDeg(li.basis, scene.true_basis, 1.0));
            std::cout << " max_angle=" << angle << " deg\n";
            expect(lili_dim == truth_dim,
                   scene.name + ": detected dimension matches truth");
            expect(angle <= 5.0,
                   scene.name + ": LiLi subspace within 5 deg of truth (H1)");
        }
    }

    // 2. Eq. (14) sparsification must keep the spanned subspace and the
    //    orthonormality of the basis, and must not increase the L1 norm.
    {
        Eigen::Matrix<double, 6, 2> U;
        U.col(0) << 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
        U.col(1) << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        const Eigen::Matrix2d mix =
            (Eigen::Matrix2d() << 0.6, -0.8, 0.8, 0.6).finished();
        U = U * mix;
        const Eigen::MatrixXd Q = U.householderQr().householderQ();
        const Eigen::MatrixXd B0 = Q.leftCols(2);
        const Eigen::Matrix<double, 6, Eigen::Dynamic> S =
            l1SparsifyBasis(U, 1.0, 50);
        const double angle = maxAngle(principalAnglesDeg(S, U, 1.0));
        expect(angle < 1e-6,
               "l1 sparsification preserves the spanned subspace");
        expect((S.transpose() * S - Eigen::Matrix2d::Identity()).norm() < 1e-9,
               "sparsified basis stays orthonormal");
        expect(S.cwiseAbs().sum() <= B0.cwiseAbs().sum() + 1e-9,
               "sparsification does not increase the L1 norm");
    }

    // 3. The perturbation set is the three axis-aligned translation+rotation
    //    couples, emitted with both signs.
    {
        PerturbationScale scale;
        scale.translation = 0.1;
        scale.rotation = 0.2;
        const std::vector<Vector6d> dirs = generatePerturbations(scale, true);
        expect(static_cast<int>(dirs.size()) == 6,
               "six bidirectional axis-coupled perturbations");
        bool coupled = true;
        for (const Vector6d& v : dirs) {
            const bool same_axis =
                (std::abs(v(0)) > 0.0 && std::abs(v(3)) > 0.0) ||
                (std::abs(v(1)) > 0.0 && std::abs(v(4)) > 0.0) ||
                (std::abs(v(2)) > 0.0 && std::abs(v(5)) > 0.0);
            coupled = coupled && same_axis;
        }
        expect(coupled,
               "perturbations couple translation and rotation on one axis");
    }

    // 4. H3: the data-reassociation mechanism is load-bearing. On a noisy
    //    non-degenerate scene the freshly-reassociated LiLi must stay quiet
    //    while the frozen-correspondence ablation hallucinates degeneracy.
    {
        SceneOptions opt;
        opt.noise_sigma = 0.02;
        const Scene random = makeRandomScene(opt);
        const IcpResult icp =
            pointToPlaneIcp(random.local, random.reference, random.P_init);
        const DetectorResult li = liliDetector(random.local, random.reference,
                                                icp.pose, icp.information);

        LiliOptions ablated;
        ablated.icp_options.reassociate = false;
        const DetectorResult li_ab = liliDetector(
            random.local, random.reference, icp.pose, icp.information, ablated);
        std::cout << "  H3 random_surface: lili_dim=" << li.basis.cols()
                  << " ablated_dim=" << li_ab.basis.cols() << "\n";
        expect(li.basis.cols() == 0,
               "H3: LiLi reports no degeneracy on a noisy random surface");
        expect(li_ab.basis.cols() > li.basis.cols(),
               "H3: freezing correspondences induces a false positive");
    }

    if (failures == 0) {
        std::cout << "LiLi detector tests passed.\n";
    }
    return failures == 0 ? 0 : 1;
}
