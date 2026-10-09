# LiLi_Test

A compact C++/Eigen reference implementation for computing a 6×6 information matrix from point-to-plane residuals and detecting pose-estimation degeneracy using eigenvalue analysis.

## Overview

This project reflects the core idea behind Lie-theory-based degeneracy detection in LiDAR scan alignment:

- each point-to-plane correspondence produces a 1×6 Jacobian in the SE(3) tangent space,
- accumulating J^T J yields the information matrix,
- eigenvalue inspection reveals near-singular directions (degenerate modes),
- the smallest eigenmodes are classified as rotation-dominated, translation-dominated, or mixed.

## Repository structure

- `include/info_matrix.h`: API and data structures.
- `src/info_matrix.cpp`: point-to-plane information accumulation and degeneracy analysis.
- `src/synthetic.cpp`: synthetic validation scenes for planar and non-planar geometry.
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

## Jacobian used

The point-to-plane residual is

r = n^T (R p + t - q)

and in the SE(3) tangent space the Jacobian is

J = [ n^T (-R [p]_x),  n^T ]

This is intentionally simple and reference-oriented: it focuses on degeneracy detection logic rather than a full SLAM pipeline.
