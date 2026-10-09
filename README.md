# LiLi_Test

A C++/Eigen reference implementation for point-to-plane residual-based information analysis and LiLi degeneracy detection.

## Overview

This repository implements the core ingredients of a LiLi-style degeneracy detector for 3D LiDAR scan alignment:

- point-to-plane residuals
- Jacobians in the SE(3) tangent space
- information matrix accumulation
- eigenvalue analysis of the 6×6 information matrix
- classification of degenerate modes as rotation-dominated, translation-dominated, or mixed

## Project structure

- `include/info_matrix.h`: API
- `src/info_matrix.cpp`: core implementation
- `src/synthetic.cpp`: synthetic planar/random-scene validation
- `tests/test_info_matrix.cpp`: smoke test
- `CMakeLists.txt`: build config

## Build

```bash
mkdir -p build
cd build
cmake ..
cmake --build .
```

## Run

```bash
./lili_synthetic
./lili_scenes                 # M1 scene acceptance report (--noise 0.02 for noise)
./lili_eval                   # M2: shared ICP + Zhang baseline vs analytic truth
ctest --output-on-failure
```

The `lili_scenes` tool checks that the analytic degeneracy basis of each
synthetic scene keeps the extended scan (Eq. 8-10) on the surface, and that a
control direction off that basis does not. `lili_eval` registers every scene
with the shared point-to-plane ICP and scores the Hessian-based (Zhang) detector
against the analytic subspace. See `docs/validation_plan.md`.

## Residual model

For a point-to-plane correspondence,

r = n^T (R p + t - q)

and the Jacobian wrt the SE(3) tangent vector is

J = [ n^T (-R [p]_x), n^T ]

Accumulating J^T J yields the information matrix used by LiLi.

## Real-world extension

This project is intentionally compact and reference-oriented. The natural next step is to integrate this into a full point-to-plane ICP pipeline and trigger a reject/downweight/update policy when the information matrix is near-singular.

## Relation to the LiLi paper

The reference paper is `docs/2609.17145v2.pdf`. This repository covers only the
static, Hessian-style part of the analysis: residual Jacobians, information-matrix
accumulation, and eigenvalue/condition-number degeneracy classification.

The paper's contribution is a *perturbation-based* detector that is not implemented
here yet: adaptive perturbation scaling from the re-association parameter `k`,
per-perturbation re-optimization, `T_degeneracy = P_opt^-1 · P_perturbed`, twist
extraction via `log(T)`, a PCA degeneracy-subspace basis with ℓ1 sparsification, and
the extended-scan alignment-quality metric `Q` (Eq. 8–10).

See `docs/validation_plan.md` for the plan to independently validate the paper's
algorithm on synthetic data.
