// M4 driver: the full experiment matrix.
//   Sweep A - noise scan over seeds, comparing Zhang / LiLi / LiLi-ablated.
//   Sweep B - perturbation-scale k scan at a fixed noise level.
//   Sweep C - threshold sensitivity (tau_PCA, tau_displacement).
// Emits one CSV row per (config, seed, scene, method, parameter) measurement.

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "detectors.h"
#include "eval_quality.h"
#include "icp.h"
#include "lili_detector.h"
#include "scenes.h"

namespace {

using lili::DetectorResult;
using lili::IcpOptions;
using lili::LiliOptions;
using lili::Scene;
using lili::SceneOptions;

struct Args {
    std::string out = "-";
    int seeds = 3;
    int q_samples = 16;
    double local_scale = 1.0;
    double local_phase = 0.0;
    double noise_for_k = 0.02;
    bool run_noise = true;
    bool run_k = true;
    bool run_thresholds = true;
};

std::string fmt(double v) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(6) << v;
    return os.str();
}

double maxAngleDeg(const Eigen::Matrix<double, 6, Eigen::Dynamic>& basis,
                   const Eigen::Matrix<double, 6, Eigen::Dynamic>& truth,
                   double length_scale) {
    if (truth.cols() == 0) { return 0.0; }
    if (basis.cols() == 0) { return 90.0; }
    const std::vector<double> angles =
        lili::principalAnglesDeg(basis, truth, length_scale);
    return angles.empty() ? 90.0 : angles.front();
}

double charLength(const lili::PointCloud& cloud) {
    if (cloud.empty()) { return 1.0; }
    std::vector<double> d;
    d.reserve(cloud.size());
    for (const auto& p : cloud.points) d.push_back(p.norm());
    std::sort(d.begin(), d.end());
    return std::max(1e-6, d[d.size() / 2]);
}

struct Recorder {
    std::ostream* out;
    lili::ExtensionOptions ext;

    void header() {
        *out << "sweep,config,seed,scene,noise,k,tau_pca,tau_disp,method,true_dim,"
                "dim,max_angle_deg,Q\n";
    }

    void emit(const std::string& sweep, const std::string& config, unsigned seed,
              const Scene& scene, double noise, int k, double tau_pca, double tau_disp,
              const std::string& method,
              const Eigen::Matrix<double, 6, Eigen::Dynamic>& basis,
              const lili::SE3& pose, double length_scale) {
        *out << sweep << ',' << config << ',' << seed << ',' << scene.name << ','
             << fmt(noise) << ',' << k << ',' << fmt(tau_pca) << ',' << fmt(tau_disp)
             << ',' << method << ',' << scene.true_basis.cols() << ',' << basis.cols()
             << ',' << fmt(maxAngleDeg(basis, scene.true_basis, length_scale)) << ','
             << fmt(lili::alignmentQuality(scene.local, scene.reference, pose, basis,
                                           ext)
                        .median)
             << '\n';
    }
};

// Runs the three methods on one scene instance and records each.
void evaluateInstance(Recorder& rec, const std::string& sweep,
                      const std::string& config, unsigned seed, const Scene& scene,
                      double noise, int k, double tau_pca, double tau_disp,
                      bool with_ablation) {
    const IcpOptions icp_opt;
    const lili::IcpResult icp =
        lili::pointToPlaneIcp(scene.local, scene.reference, scene.P_init, icp_opt);
    const double L = charLength(scene.local);

    const DetectorResult zhang = lili::zhangDetector(icp.information);
    rec.emit(sweep, config, seed, scene, noise, k, tau_pca, tau_disp, "zhang",
             zhang.basis, icp.pose, L);

    LiliOptions lo;
    lo.reassociation_k = k;
    lo.pca_rel_threshold = tau_pca;
    lo.tau_displacement = tau_disp;
    const DetectorResult li =
        lili::liliDetector(scene.local, scene.reference, icp.pose, icp.information, lo);
    rec.emit(sweep, config, seed, scene, noise, k, tau_pca, tau_disp, "lili",
             li.basis, icp.pose, L);

    if (with_ablation) {
        LiliOptions ab = lo;
        ab.icp_options.reassociate = false;
        const DetectorResult li_ab = lili::liliDetector(
            scene.local, scene.reference, icp.pose, icp.information, ab);
        rec.emit(sweep, config, seed, scene, noise, k, tau_pca, tau_disp,
                 "lili_ablated", li_ab.basis, icp.pose, L);
    }
}

Args parse(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--out" && i + 1 < argc) a.out = argv[++i];
        else if (arg == "--seeds" && i + 1 < argc) a.seeds = std::atoi(argv[++i]);
        else if (arg == "--q-samples" && i + 1 < argc) a.q_samples = std::atoi(argv[++i]);
        else if (arg == "--local-scale" && i + 1 < argc) a.local_scale = std::atof(argv[++i]);
        else if (arg == "--local-phase" && i + 1 < argc) a.local_phase = std::atof(argv[++i]);
        else if (arg == "--noise" && i + 1 < argc) a.noise_for_k = std::atof(argv[++i]);
        else if (arg == "--only-noise") { a.run_k = false; a.run_thresholds = false; }
        else if (arg == "--only-k") { a.run_noise = false; a.run_thresholds = false; }
        else if (arg == "--only-thresholds") { a.run_noise = false; a.run_k = false; }
        else {
            std::cerr << "Usage: " << argv[0]
                      << " [--out F] [--seeds n] [--q-samples n] [--local-scale s]"
                      << " [--local-phase p] [--noise s]"
                      << " [--only-noise|--only-k|--only-thresholds]\n";
            std::exit(2);
        }
    }
    return a;
}

}  // namespace

int main(int argc, char** argv) {
    const Args args = parse(argc, argv);

    std::ofstream file;
    std::ostream* out = &std::cout;
    if (args.out != "-") {
        file.open(args.out);
        if (!file) {
            std::cerr << "cannot open " << args.out << "\n";
            return 1;
        }
        out = &file;
    }

    Recorder rec;
    rec.out = out;
    rec.ext.n_samples = args.q_samples;
    rec.header();

    const std::string config =
        "scale" + fmt(args.local_scale) + "_phase" + fmt(args.local_phase);

    if (args.run_noise) {
        const std::vector<double> noise_levels = {0.0, 0.01, 0.02, 0.03, 0.05};
        for (double noise : noise_levels) {
            for (int s = 0; s < args.seeds; ++s) {
                SceneOptions opt;
                opt.noise_sigma = noise;
                opt.local_spacing_scale = args.local_scale;
                opt.local_phase = args.local_phase;
                opt.seed = 20261009u + 1013u * static_cast<unsigned>(s);
                const unsigned seed = opt.seed;
                for (const Scene& scene : lili::makeAllScenes(opt)) {
                    evaluateInstance(rec, "noise", config, seed, scene, noise,
                                     /*k=*/3, /*tau_pca=*/3e-2, /*tau_disp=*/0.10,
                                     /*with_ablation=*/true);
                }
            }
        }
    }

    if (args.run_k) {
        const std::vector<int> ks = {1, 3, 5, 10};
        for (int k : ks) {
            for (int s = 0; s < args.seeds; ++s) {
                SceneOptions opt;
                opt.noise_sigma = args.noise_for_k;
                opt.local_spacing_scale = args.local_scale;
                opt.local_phase = args.local_phase;
                opt.seed = 20261009u + 1013u * static_cast<unsigned>(s);
                const unsigned seed = opt.seed;
                for (const Scene& scene : lili::makeAllScenes(opt)) {
                    evaluateInstance(rec, "k", config, seed, scene, args.noise_for_k, k,
                                     3e-2, 0.10, /*with_ablation=*/false);
                }
            }
        }
    }

    if (args.run_thresholds) {
        const std::vector<double> pca_rels = {1e-2, 3e-2, 1e-1};
        const std::vector<double> taus = {0.05, 0.10, 0.30};
        for (double pr : pca_rels) {
            for (double td : taus) {
                SceneOptions opt;
                opt.noise_sigma = 0.02;
                opt.local_spacing_scale = args.local_scale;
                opt.local_phase = args.local_phase;
                opt.seed = 20261009u;
                for (const Scene& scene : lili::makeAllScenes(opt)) {
                    evaluateInstance(rec, "threshold", config, opt.seed, scene, 0.02,
                                     /*k=*/3, pr, td, /*with_ablation=*/false);
                }
            }
        }
    }

    if (out == &std::cout) {
        // Nothing else to do; CSV already on stdout.
    } else {
        std::cerr << "wrote " << args.out << "\n";
    }
    return 0;
}
