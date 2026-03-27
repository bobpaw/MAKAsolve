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

Ensure that PUMI is installed somewhere that is in the `CMAKE_PREFIX_PATH`.
Assuming PUMI is installed in `/opt/PUMI`, the following example `build.sh`
can be run from the project folder to configure and build MAKAsolve:

```sh
CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/opt/PUMI
cmake -S . -B build
cmake --build build
```

## Authors

- Mikiel Gica <gicam@rpi.edu>
- Aidan Hoover <hoovea@rpi.edu>
- Kiran Koch-Mathews <kochmk@rpi.edu>
- Aiden Woodruff <woodra@rpi.edu>

[CMake]: https://cmake.org
[PUMI]: https://github.com/SCOREC/core
