# adFE2D
A finite element model of the advection diffusion equation

## Installation
To compile adFE2D, just use `make` in the project directory. So far
has only been tested on Linux, and probably won't work on Windows
because it uses the `shell pwd` command  to determine the project
directory.

Required Programs:
* `gfortran` to compile the fortran source code
* `gmsh` to create input mesh file
* `BLAS` and `LAPACK` for linear algebra

Optional python3 modules for visualization:
* STILL NEEDS TO FILLED OUT

Plus a sample image:
[SAMPLE IMAGE]
