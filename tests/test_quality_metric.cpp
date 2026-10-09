#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "eval_quality.h"
#include "scenes.h"

namespace {

int failures = 0;

void expect(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    }
}

// One se(3) direction outside the true degeneracy subspace.
Eigen::Matrix<double, 6, Eigen::Dynamic> controlDirection(
    const Eigen::Matrix<double, 6, Eigen::Dynamic>& truth) {
    const Eigen::Matrix<double, 6, 6> Q = truth.householderQr().householderQ();
    return Q.rightCols(6 - static_cast<int>(truth.cols())).col(0);
}

}  // namespace

int main() {
    using namespace lili;

    // 1. A cloud compared against itself has zero distance everywhere.
    {
        SceneOptions opt;
        opt.noise_sigma = 0.0;
        const Scene scene = makePlaneScene(opt);
        const auto q = alignmentQualitySingle(scene.local, scene.local, SE3::identity());
        expect(q.median < 1e-12, "identical clouds give Q == 0");
    }

    // 2. With the exact pose and an empty basis, Q is at the noise floor.
    // 3. The analytic degeneracy basis keeps Q near that floor, while a control
    //    direction outside the basis degrades it substantially.
    {
        SceneOptions opt;
        opt.noise_sigma = 0.0;
        for (const Scene& scene : makeAllScenes(opt)) {
            const auto q_plain =
                alignmentQualitySingle(scene.local, scene.reference, scene.P_gt);
            const auto q_true =
                alignmentQuality(scene.local, scene.reference, scene.P_gt, scene.true_basis);
            const auto q_control = alignmentQuality(
                scene.local, scene.reference, scene.P_gt, controlDirection(scene.true_basis));

            std::cout << scene.name << ": Q_plain=" << q_plain.median
                      << " Q_true=" << q_true.median
                      << " Q_control=" << q_control.median << "\n";

            expect(q_plain.median < 1e-6,
                   scene.name + ": exact pose aligns with Q ~ 0");
            expect(q_true.median < q_control.median,
                   scene.name + ": true basis beats the control direction");
        }
    }

    // 4. Principal angles: a basis against itself is 0 degrees, and orthogonal
    //    subspaces are 90 degrees apart.
    {
        Eigen::Matrix<double, 6, Eigen::Dynamic> A(6, 2);
        A.setZero();
        A(0, 0) = 1.0;  // rotation x
        A(1, 1) = 1.0;  // rotation y
        const auto same = principalAnglesDeg(A, A);
        expect(!same.empty() && same.front() < 1e-6, "principal angle to itself is 0");

        Eigen::Matrix<double, 6, Eigen::Dynamic> B(6, 2);
        B.setZero();
        B(2, 0) = 1.0;  // rotation z
        B(3, 1) = 1.0;  // translation x
        const auto ort = principalAnglesDeg(A, B);
        expect(!ort.empty() && std::abs(ort.front() - 90.0) < 1e-6,
               "orthogonal subspaces are 90 degrees apart");
    }

    if (failures == 0) {
        std::cout << "LiLi quality-metric tests passed.\n";
    }
    return failures == 0 ? 0 : 1;
}
