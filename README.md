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
ctest --output-on-failure
```

## Residual model

For a point-to-plane correspondence,

r = n^T (R p + t - q)

and the Jacobian wrt the SE(3) tangent vector is

J = [ n^T (-R [p]_x), n^T ]

Accumulating J^T J yields the information matrix used by LiLi.

## Real-world extension

This project is intentionally compact and reference-oriented. The natural next step is to integrate this into a full point-to-plane ICP pipeline and trigger a reject/downweight/update policy when the information matrix is near-singular.
