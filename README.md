# MAKAsolve

MAKAsolve is a Parallel GPU (CUDA) Finite Element Stabilized
Advection-Diffusion solver.

## Requirements

- CMake >= 3.12
- [PUMI](github.com/SCOREC/core) aka SCOREC core >= 4
- [HYPRE](https://github.com/hypre-space/hypre) >= 3.1.0
- MPI
- CUDA >= 11.2

## Build instructions

1. Install dependencies
2. Clone MAKAsolve
2. Configure CMake
3. Build

Ensure that PUMI is installed somewhere that is in the `CMAKE_PREFIX_PATH`.
Assuming PUMI and HYPRE are installed in `/opt/PUMI` and `/opt/HYPRE`,
the following example `build.sh` can be run from the project folder to
configure and build MAKAsolve:

```sh
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/opt/PUMI:/opt/HYPRE
cmake -S . -B build
cmake --build build
```

### Installing HYPRE on AiMOS

HYPRE can be built with GPU support on AiMOS with the 
[NVIDIA HPC SDK](https://docs.cci.rpi.edu/software/NVHPC/) for CUDA >= 11.2. 
HYPRE can be built and installed on AiMOS with CMake 
[(link to instructions)](https://github.com/hypre-space/hypre/blob/master/INSTALL.md#hypre-installation-information-using-cmake) if HYPRE's minimum CMake version 
is lowered (this does not appear to cause any issues). `HYPRE_ENABLE_CUDA` must be enabled.

## Authors

- Mikiel Gica <gicam@rpi.edu>
- Aidan Hoover <hoovea@rpi.edu>
- Kiran Koch-Mathews <kochmk@rpi.edu>
- Aiden Woodruff <woodra@rpi.edu>

[CMake]: https://cmake.org
[PUMI]: https://github.com/SCOREC/core
