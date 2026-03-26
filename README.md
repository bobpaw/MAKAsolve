# MAKAsolve

MAKAsolve is a Parallel GPU (CUDA) Finite Element Stabilized
Advection-Diffusion solver.

## Requirements

- [CMake]() >= 3.12
- [PUMI]() aka SCOREC core >= 4
- MPI
- CUDA

## Build instructions

1. Install dependencies
2. Clone MAKAsolve
2. Configure CMake
3. Build

An example `build.sh` that can be run from the project folder is given below:

```sh
CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/opt/PUMI
cmake -S . -B build
cmake --build build
```

## Authors

- Aiden Woodruff <woodra@rpi.edu>
- Mikiel Gica <gicam@rpi.edu>
- Kiran Koch-Mathews <kochmk@rpi.edu>
- Aidan Hoover <hoovea@rpi.edu>

[CMake]: https://cmake.org
[PUMI]: https://github.com/SCOREC/core
