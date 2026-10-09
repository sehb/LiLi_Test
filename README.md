# LiLi_Test

A compact C++/Eigen reference implementation for computing a 6x6 information matrix from point-to-plane residuals and detecting pose-estimation degeneracy using eigenvalue analysis.

## Overview

This project is designed to reflect the core idea behind Lie-theory-based degeneracy detection in LiDAR scan alignment:

- each point-to-plane correspondence produces a 1x6 Jacobian in the SE(3) tangent space,
- accumulating J^T J yields the information matrix,
- eigenvalue inspection reveals near-singular directions (degenerate modes),
- the smallest eigenmodes are classified as rotation-dominated, translation-dominated, or mixed.

## Repository structure

- `include/info_matrix.h`: API and data structures.
- `src/info_matrix.cpp`: point-to-plane information accumulation and degeneracy analysis.
- `src/synthetic.cpp`: synthetic validation scenes for planar and non-planar motions.
- `tests/test_info_matrix.cpp`: smoke tests for the information matrix.
- `CMakeLists.txt`: build configuration.

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
ctest --output-on-failure
```

## Notes

The Jacobian used here is based on the residual

r = n^T (R p + t - q)

with the SE(3) perturbation expressed in the tangent space, giving

J = [ n^T (-R [p]_x),  n^T ]

This is intentionally simple and reference-oriented: it focuses on the degeneracy detection logic rather than a full SLAM pipeline.
