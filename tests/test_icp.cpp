#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "detectors.h"
#include "eval_quality.h"
#include "icp.h"
#include "scenes.h"

namespace {

int failures = 0;

void expect(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << "\n";
        ++failures;
    }
}

double rotationDeg(const lili::SE3& a, const lili::SE3& b) {
    const lili::SE3 e = a.inverse() * b;
    const double c = std::max(-1.0, std::min(1.0, (e.R.trace() - 1.0) * 0.5));
    return std::acos(c) * 180.0 / M_PI;
}

// Characteristic length used to make rotation and translation comparable.
double charLength(const lili::Scene& scene) {
    std::vector<double> d;
    d.reserve(scene.local.size());
    for (const auto& p : scene.local.points) { d.push_back(p.norm()); }
    if (d.empty()) { return 1.0; }
    std::sort(d.begin(), d.end());
    return std::max(1e-6, d[d.size() / 2]);
}

}  // namespace

int main() {
    using namespace lili;

    // 1. An identity copy must converge immediately with (near) zero residual.
    {
        SceneOptions opt;
        const Scene scene = makePlaneScene(opt);
        const IcpResult r = pointToPlaneIcp(scene.local, scene.reference, scene.P_gt);
        expect(r.converged, "identity start converges");
        expect(r.rms < 1e-9, "identity start has zero residual");
        expect(r.correspondences == static_cast<int>(scene.local.size()),
               "identity start matches every local point");
    }

    // 2. Registration on a degenerate scene is not unique: the correct
    //    noise-free statement is that the difference between the converged pose
    //    and the ground truth lies inside the known degeneracy subspace. On the
    //    non-degenerate control the pose must be recovered outright.
    {
        SceneOptions opt;
        opt.noise_sigma = 0.0;
        for (const Scene& scene : makeAllScenes(opt)) {
            const IcpResult r = pointToPlaneIcp(scene.local, scene.reference, scene.P_init);
            const double rot = rotationDeg(r.pose, scene.P_gt);
            const double trans = (r.pose.t - scene.P_gt.t).norm();
            const Vector6d err = logSE3(scene.P_gt.inverse() * r.pose);
            const auto dev = deviationFromSubspace(err, scene.true_basis, charLength(scene));
            std::cout << "  " << scene.name << ": iter=" << r.iterations
                      << " rot_err=" << rot << " trans_err=" << trans
                      << " rms=" << r.rms << " |err|=" << dev.total
                      << " orthogonal=" << dev.orthogonal
                      << " frac=" << dev.fraction << "\n";
            expect(r.converged, scene.name + ": ICP converged");
            if (scene.true_basis.cols() == 0) {
                // Fully constrained: nothing may absorb the error.
                expect(dev.total < 0.02,
                       scene.name + ": non-degenerate scene recovers the pose");
            } else {
                // The ambiguity must be explained by the known subspace.
                expect(dev.fraction < 0.10,
                       scene.name + ": pose error lies in the true degeneracy subspace");
            }
        }
    }

    // 2b. Under noise the pose error is no longer confined to the analytic
    //     subspace: noise excites every weakly constrained direction. The
    //     checkable property is that registration stays bounded and converges.
    {
        SceneOptions opt;
        opt.noise_sigma = 0.02;
        for (const Scene& scene : makeAllScenes(opt)) {
            const IcpResult r = pointToPlaneIcp(scene.local, scene.reference, scene.P_init);
            const Vector6d err = logSE3(scene.P_gt.inverse() * r.pose);
            const double total =
                deviationFromSubspace(err, scene.true_basis, charLength(scene)).total;
            std::cout << "  " << scene.name << " (noisy): it=" << r.iterations
                      << " rms=" << r.rms << " |err|=" << total << "\n";
            expect(r.converged, scene.name + ": noisy ICP converged");
            expect(total < 0.25, scene.name + ": noisy pose error stays bounded");
        }
    }

    // 3. The information matrix returned by the ICP must be symmetric PSD, and
    //    the detector must see the plane's degeneracy but not a plane's worth of
    //    it on the non-degenerate control.
    {
        SceneOptions opt;
        const Scene plane = makePlaneScene(opt);
        const IcpResult rp = pointToPlaneIcp(plane.local, plane.reference, plane.P_init);
        expect((rp.information - rp.information.transpose()).norm() < 1e-9,
               "information matrix is symmetric");
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> solver(
            rp.information);
        expect(solver.eigenvalues().minCoeff() > -1e-9,
               "information matrix is positive semidefinite");

        const DetectorResult zhang = zhangDetector(rp.information);
        expect(static_cast<int>(plane.true_basis.cols()) ==
                   static_cast<int>(zhang.basis.cols()),
               "Zhang finds 3 degenerate directions on the plane");

        const Scene random = makeRandomScene(opt);
        const IcpResult rr =
            pointToPlaneIcp(random.local, random.reference, random.P_init);
        const DetectorResult zhang_random = zhangDetector(rr.information);
        expect(zhang_random.basis.cols() == 0,
               "Zhang reports no degeneracy on the non-degenerate control");
    }

    if (failures == 0) {
        std::cout << "LiLi ICP / Zhang tests passed.\n";
    }
    return failures == 0 ? 0 : 1;
}
