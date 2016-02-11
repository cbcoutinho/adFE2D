# adFE2D
A finite element model of the advection diffusion equation

## Installation
To compile adFE2D, just use `make` in the project directory. So far
the makefile has only been tested on Linux, and probably won't work
on Windows because it uses the `shell pwd` command  to determine the
project directory.

Required Programs:
* `gfortran` to compile the fortran source code
* `gmsh` to create input mesh file
* `BLAS` and `LAPACK` for linear algebra

Optional python3 modules for visualization:
* STILL NEEDS TO FILLED OUT

Plus a sample image:
[SAMPLE IMAGE]

## To run
To run adFE2D, you first need to create a mesh file using GMSH, or
manually write a .msh file that has the same format as a GMSH .msh
file. Currently, only linear quadrilateral meshes are allowed, but
linear trianle element meshes as well as hybrid meshes are being
implemented.

Once the the program is compiled and a mesh file has been created,
the program is executed using the following command:

`../bin/adFE2D example.msh`

The output of adFE2D will be sent to:

* 'stdout' with time and convergence information
* 'data.dat' with xy-cooridnates and output data in columns (based on frequency in params.in)
* 'scratch.dat' with convergence data

## Plot residuals from scratch.dat using gnuplot:
`gnuplot -persist -e "plot 'scratch.dat'" loop.plt`
