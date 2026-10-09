#include "scenes.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace lili {
namespace {

struct Sample {
    Eigen::Vector3d p;
    Eigen::Vector3d n;
};

struct Shape {
    std::vector<Sample> samples;
    void push(const Sample& s) { samples.push_back(s); }
    std::size_t size() const { return samples.size(); }
};

// Reference and local view of one surface. They are separate samplings: a real
// scan is not a subset of the map it is registered against.
struct Surface {
    Shape reference;
    Shape local;
};

constexpr double kPi = 3.14159265358979323846;

void thinShape(Shape& shape, int max_points) {
    if (max_points <= 0 || shape.size() <= static_cast<std::size_t>(max_points)) {
        return;
    }
    const std::size_t stride = shape.size() / static_cast<std::size_t>(max_points);
    if (stride < 2) { return; }
    Shape thinned;
    thinned.samples.reserve(shape.size() / stride + 1);
    for (std::size_t i = 0; i < shape.size(); i += stride) {
        thinned.push(shape.samples[i]);
    }
    shape = thinned;
}

// Builds the two samplings of a surface parametrised by (u, v).
//
// With the default options the local view is selected from the *same* grid as
// the reference, so P_gt maps every local point exactly onto a reference point
// and the noise-free alignment floor is zero. Changing `local_spacing_scale`
// and/or `local_phase` instead samples the local window on its own grid, which
// is what makes nearest-neighbour association ambiguous.
template <typename EmitFn>
Surface gridSurface(double u_lo, double u_hi, double u_step, double v_lo,
                    double v_hi, double v_step, double u_loc_lo, double u_loc_hi,
                    double v_loc_lo, double v_loc_hi, EmitFn emit,
                    const SceneOptions& opt) {
    Surface surface;
    std::vector<char> is_local;

    for (double u = u_lo; u <= u_hi + 1e-9; u += u_step) {
        for (double v = v_lo; v <= v_hi + 1e-9; v += v_step) {
            surface.reference.push(emit(u, v));
            is_local.push_back((u >= u_loc_lo - 1e-9 && u <= u_loc_hi + 1e-9 &&
                                v >= v_loc_lo - 1e-9 && v <= v_loc_hi + 1e-9)
                                   ? 1
                                   : 0);
        }
    }

    const bool independent =
        opt.local_spacing_scale != 1.0 || opt.local_phase != 0.0;
    if (!independent) {
        for (std::size_t i = 0; i < surface.reference.size(); ++i) {
            if (is_local[i]) { surface.local.push(surface.reference.samples[i]); }
        }
        return surface;
    }

    const double u_step_local = u_step * opt.local_spacing_scale;
    const double v_step_local = v_step * opt.local_spacing_scale;
    const double u_start = u_loc_lo + opt.local_phase * u_step_local;
    const double v_start = v_loc_lo + opt.local_phase * v_step_local;
    for (double u = u_start; u <= u_loc_hi + 1e-9; u += u_step_local) {
        for (double v = v_start; v <= v_loc_hi + 1e-9; v += v_step_local) {
            surface.local.push(emit(u, v));
        }
    }
    return surface;
}

void addNoise(PointCloud& cloud, double sigma, std::mt19937_64& rng) {
    if (sigma <= 0.0) { return; }
    std::normal_distribution<double> dist(0.0, sigma);
    for (auto& p : cloud.points) {
        p += Eigen::Vector3d(dist(rng), dist(rng), dist(rng));
    }
}

PointCloud toWorld(const Shape& shape) {
    const Eigen::Matrix3d Rw = shapeToWorldRotation();
    PointCloud cloud;
    cloud.reserve(shape.size());
    for (const auto& s : shape.samples) {
        cloud.points.push_back(Rw * s.p);
        cloud.normals.push_back(Rw * s.n);
    }
    return cloud;
}

PointCloud toSensorFrame(const Shape& shape, const Eigen::Vector3d& sensor) {
    PointCloud cloud;
    cloud.reserve(shape.size());
    for (const auto& s : shape.samples) {
        cloud.points.push_back(s.p - sensor);
        cloud.normals.push_back(s.n);
    }
    return cloud;
}

Eigen::Matrix<double, 6, Eigen::Dynamic> makeBasis(
    const std::vector<Vector6d>& twists) {
    Eigen::Matrix<double, 6, Eigen::Dynamic> B(6, static_cast<int>(twists.size()));
    for (std::size_t i = 0; i < twists.size(); ++i) {
        B.col(static_cast<int>(i)) = twists[i];
    }
    return B;
}

Scene assemble(const std::string& name, Surface surface,
               const Eigen::Vector3d& sensor_shape,
               const std::vector<Vector6d>& shape_twists, bool confirmed,
               const std::string& notes, const std::string& expected,
               const SceneOptions& opt) {
    Scene scene;
    scene.name = name;
    scene.notes = notes;
    scene.expected_comment = expected;
    scene.ground_truth_confirmed = confirmed;

    thinShape(surface.reference, opt.max_points);
    thinShape(surface.local, opt.max_points);

    const Eigen::Matrix3d Rw = shapeToWorldRotation();
    scene.reference = toWorld(surface.reference);
    scene.local = toSensorFrame(surface.local, sensor_shape);
    scene.P_gt = SE3{Rw, Rw * sensor_shape};

    std::vector<Vector6d> body_twists;
    body_twists.reserve(shape_twists.size());
    for (const auto& xi : shape_twists) {
        body_twists.push_back(twistAtPoint(xi, sensor_shape));
    }
    scene.true_basis = makeBasis(body_twists);

    Vector6d delta0;
    delta0 << 0.015, -0.010, 0.012, 0.030, -0.020, 0.025;
    scene.P_init = opt.apply_initial_offset ? expSE3(delta0) * scene.P_gt : scene.P_gt;

    std::mt19937_64 rng(opt.seed);
    addNoise(scene.reference, opt.noise_sigma, rng);
    addNoise(scene.local, opt.noise_sigma, rng);
    return scene;
}

// --- shape emitters, all parametrised as emit(u, v) -> Sample ---------------

Sample planeEmit(double x, double y) {
    return Sample{Eigen::Vector3d(x, y, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0)};
}

struct CylinderEmit {
    double radius;
    Sample operator()(double z, double th) const {
        const Eigen::Vector3d outward(std::cos(th), std::sin(th), 0.0);
        return Sample{radius * outward + Eigen::Vector3d(0, 0, z), -outward};
    }
};

struct SinusoidalEmit {
    double r0, amp, omega;
    Sample operator()(double z, double th) const {
        const double r = r0 + amp * std::sin(omega * z);
        const double dr = amp * omega * std::cos(omega * z);
        const Eigen::Vector3d p(r * std::cos(th), r * std::sin(th), z);
        const Eigen::Vector3d n =
            Eigen::Vector3d(std::cos(th), std::sin(th), -dr).normalized();
        return Sample{p, -n};
    }
};

double roughHeight(double x, double y) {
    return 0.14 * std::sin(3.1 * x + 0.7) * std::cos(2.3 * y - 1.1) +
           0.10 * std::sin(5.2 * y + 0.3) +
           0.08 * std::cos(4.1 * x - 0.5) +
           0.06 * std::sin(6.7 * x + 1.3) * std::sin(7.9 * y - 0.9);
}

void roughGradient(double x, double y, double* dzdx, double* dzdy) {
    const double s1 = std::sin(3.1 * x + 0.7), c1 = std::cos(3.1 * x + 0.7);
    const double s2 = std::cos(2.3 * y - 1.1), c2 = -std::sin(2.3 * y - 1.1);
    const double s4 = std::sin(6.7 * x + 1.3), c4 = std::cos(6.7 * x + 1.3);
    const double s5 = std::sin(7.9 * y - 0.9), c5 = std::cos(7.9 * y - 0.9);
    *dzdx = 0.14 * 3.1 * c1 * s2 - 0.08 * 4.1 * std::sin(4.1 * x - 0.5) +
            0.06 * 6.7 * c4 * s5;
    *dzdy = 0.14 * 2.3 * s1 * c2 + 0.10 * 5.2 * std::cos(5.2 * y + 0.3) +
            0.06 * 7.9 * s4 * c5;
}

Sample roughEmit(double x, double y) {
    double dzdx = 0.0, dzdy = 0.0;
    roughGradient(x, y, &dzdx, &dzdy);
    return Sample{Eigen::Vector3d(x, y, roughHeight(x, y)),
                  Eigen::Vector3d(-dzdx, -dzdy, 1.0).normalized()};
}

// Screw-symmetric tube: sweep a curved (elliptical) profile with a screw
// motion, p(s, t) = (rho(s) cos(a t), rho(s) sin(a t), B sin(s) + t).
// Exactly invariant under the screw that translates along z while rotating
// about it, and - because the profile is curved, not straight - not a ruled
// surface, so there is no additional near-free sliding along rulings.
// Stands in for the paper's "Rotating Tunnel" coupled-degeneracy scenario.
struct ScrewTubeEmit {
    double R, A, B, alpha;
    Sample operator()(double s, double t) const {
        const auto point = [&](double s_, double t_) {
            const double rho = R + A * std::cos(s_);
            const double phi = alpha * t_;
            return Eigen::Vector3d(rho * std::cos(phi), rho * std::sin(phi),
                                   B * std::sin(s_) + t_);
        };
        const double h = 1e-5;
        const Eigen::Vector3d ds = (point(s + h, t) - point(s - h, t)) / (2.0 * h);
        const Eigen::Vector3d dt = (point(s, t + h) - point(s, t - h)) / (2.0 * h);
        Eigen::Vector3d n = ds.cross(dt);
        const double norm = n.norm();
        if (norm < 1e-12) { n = Eigen::Vector3d(0.0, 0.0, 1.0); }
        else { n /= norm; }
        if (n.dot(-point(s, t)) < 0.0) { n = -n; }  // face the axis
        return Sample{point(s, t), n};
    }
};

}  // namespace

Eigen::Matrix3d shapeToWorldRotation() {
    return (Eigen::AngleAxisd(0.35, Eigen::Vector3d::UnitX()) *
            Eigen::AngleAxisd(-0.20, Eigen::Vector3d::UnitY()) *
            Eigen::AngleAxisd(0.50, Eigen::Vector3d::UnitZ()))
        .toRotationMatrix();
}

Scene makePlaneScene(const SceneOptions& opt) {
    const Eigen::Vector3d sensor(0.0, 0.0, 1.0);
    Surface surface = gridSurface(-1.0, 1.0, opt.spacing, -1.0, 1.0, opt.spacing,
                                  -0.6, 0.6, -0.6, 0.6, planeEmit, opt);
    std::vector<Vector6d> twists(3, Vector6d::Zero());
    twists[0].tail<3>() = Eigen::Vector3d::UnitX();
    twists[1].tail<3>() = Eigen::Vector3d::UnitY();
    twists[2].head<3>() = Eigen::Vector3d::UnitZ();
    return assemble("plane", surface, sensor, twists, true,
                    "Bounded plane patch; boundary effects make rotation about the "
                    "normal only asymptotically degenerate.",
                    "2 in-plane translations + rotation about the normal", opt);
}

Scene makeOpenCylinderScene(const SceneOptions& opt) {
    const double radius = 1.0;
    const double d_theta = opt.spacing / radius;
    Surface surface = gridSurface(-0.6, 0.6, opt.spacing, -kPi, kPi, d_theta,
                                  -0.35, 0.35, -0.7, 0.7, CylinderEmit{radius}, opt);
    const Eigen::Vector3d sensor(0.55, 0.0, 0.0);
    std::vector<Vector6d> twists(2, Vector6d::Zero());
    twists[0].tail<3>() = Eigen::Vector3d::UnitZ();  // along the axis
    twists[1].head<3>() = Eigen::Vector3d::UnitZ();  // about the axis
    return assemble("open_cylinder", surface, sensor, twists, true,
                    "Side surface only. Axis rotation becomes a coupled "
                    "translation+rotation in the sensor frame because the sensor "
                    "is off-axis.",
                    "axis translation + axis rotation (coupled in body frame)", opt);
}

Scene makeClosedCylinderScene(const SceneOptions& opt) {
    const double radius = 1.0;
    const double d_theta = opt.spacing / radius;
    Surface surface = gridSurface(-0.6, 0.6, opt.spacing, -kPi, kPi, d_theta,
                                  -0.35, 0.35, -0.7, 0.7, CylinderEmit{radius}, opt);
    // Caps are reference-only: they never fall inside the local side window.
    const Eigen::Vector3d cap_normal(0.0, 0.0, 1.0);
    for (double z : {0.6, -0.6}) {
        for (double x = -radius; x <= radius + 1e-9; x += opt.spacing) {
            for (double y = -radius; y <= radius + 1e-9; y += opt.spacing) {
                if (x * x + y * y <= radius * radius) {
                    surface.reference.push(
                        Sample{Eigen::Vector3d(x, y, z), cap_normal});
                }
            }
        }
    }
    const Eigen::Vector3d sensor(0.55, 0.0, 0.0);
    std::vector<Vector6d> twists(2, Vector6d::Zero());
    twists[0].tail<3>() = Eigen::Vector3d::UnitZ();
    twists[1].head<3>() = Eigen::Vector3d::UnitZ();
    return assemble("closed_cylinder", surface, sensor, twists, true,
                    "Capped cylinder; the caps add no degeneracy.",
                    "axis translation + axis rotation", opt);
}

Scene makeSinusoidalCylinderScene(const SceneOptions& opt) {
    const double r0 = 1.0;
    const double d_theta = opt.spacing / r0;
    Surface surface = gridSurface(-0.6, 0.6, opt.spacing, -kPi, kPi, d_theta,
                                  -0.35, 0.35, -0.7, 0.7,
                                  SinusoidalEmit{1.0, 0.12, 2.5}, opt);
    const Eigen::Vector3d sensor(0.55, 0.0, 0.0);
    // r(z) is independent of theta, so rotation about the axis is exact while
    // axial translation is not (it maps r(z) to r(z + dz)).
    std::vector<Vector6d> twists(1, Vector6d::Zero());
    twists[0].head<3>() = Eigen::Vector3d::UnitZ();
    return assemble("sinusoidal_cylinder", surface, sensor, twists, true,
                    "Radius varies with z; only the axial rotation survives as an "
                    "exact symmetry.",
                    "axis rotation only (1-D)", opt);
}

Scene makeRandomScene(const SceneOptions& opt) {
    const Eigen::Vector3d sensor(0.0, 0.0, 1.2);
    Surface surface = gridSurface(-1.0, 1.0, opt.spacing, -1.0, 1.0, opt.spacing,
                                  -0.6, 0.6, -0.6, 0.6, roughEmit, opt);
    return assemble("random_surface", surface, sensor, {}, true,
                    "Non-degenerate control: rough surface with no invariant motion. "
                    "Reporting degeneracy here is a false positive.",
                    "none (non-degenerate control)", opt);
}

Scene makeScrewTubeScene(const SceneOptions& opt) {
    const double R = 1.2, A = 0.35, B = 0.5, alpha = 0.30;
    Surface surface = gridSurface(0.0, 2.0 * kPi, opt.spacing, -1.5, 1.5, opt.spacing,
                                  kPi - 0.8, kPi + 0.8, -0.5, 0.5,
                                  ScrewTubeEmit{R, A, B, alpha}, opt);
    // The sensor sits on the screw axis, so the body-frame twist is exactly
    // (0, 0, alpha, 0, 0, 1): a coupled rotation + translation.
    const Eigen::Vector3d sensor(0.0, 0.0, 0.0);
    Vector6d screw = Vector6d::Zero();
    screw.head<3>() = alpha * Eigen::Vector3d::UnitZ();
    screw.tail<3>() = Eigen::Vector3d::UnitZ();
    return assemble("screw_tube", surface, sensor, {screw}, true,
                    "Screw-symmetric tube with a curved (non-ruled) profile; the "
                    "exact symmetry is the screw motion itself. Stands in for the "
                    "paper's Rotating Tunnel scenario.",
                    "screw motion along the axis (1-D, coupled)", opt);
}

std::vector<Scene> makeAllScenes(const SceneOptions& opt) {
    return {makePlaneScene(opt), makeClosedCylinderScene(opt),
            makeOpenCylinderScene(opt), makeSinusoidalCylinderScene(opt),
            makeScrewTubeScene(opt), makeRandomScene(opt)};
}

}  // namespace lili
