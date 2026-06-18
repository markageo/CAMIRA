
<!-- # <img src="images/logo.png?raw=true" width="128" valign="middle" alt="RAJA"/> -->

# CAMIRA
CAMIRA (**C**oupled **A**lgorithm for **M**ultigrid **I**mmersed Boundary **R**eynolds **A**veraged Navier-Stokes, also an Aboriginal name for the word "windy") is a highly efficient solver for the incompressible steady-state Reynolds-Averaged Navier-Stokes equations. It also contains solver for calculating the unsteady concentration of a passive scalar in the flow field. The code and the methods used were developed as part of my PhD. 

CAMIRA contains two applications: 
* CAMIRA-FLOW, which is a fully coupled matrix free smoother with Full Approximation Scheme (FAS) multigrid to solve the equations on a rectillinear collocated finite volume grid. Complex geometries are accounted for using a mass conservative directional immersed boundary method. Details on these numerical methods can be found in the various [publications](#publications) that are part of this work.
* CAMIRA-PLUME, which uses Lagrangian Particle dispersion to calculate the concentration of a passive scalar in the resulting flow field given by CAMIRA-FLOW.

CAMIRA can run in parallel on shared memory CPUs using OpenMP.


## Getting Started

### Prerequisites

CAMIRA makes use of the [Eigen C++ linear algebra library](https://eigen.tuxfamily.org/index.php?title=Main_Page) and the [Computational Geometry Algorithms Library (CGAL)](https://www.cgal.org/). These are included within the source code and do not need to be installed by the user. 

The CAMIRA code uses a number of external dependencies to work. These need to be installed on the system in order to successfully build and run the code. 

* [VTK](https://vtk.org/) (Optional) - This is required to read in initial conditions for VTK files. It is not necessary to build and use the code, but required if you wish to read an initial condition field from file.
* [GMP](https://gmplib.org/), [MPFR](https://www.mpfr.org/), and [Boost](https://www.boost.org/) -  All required for [CGAL](https://www.cgal.org/).
* [CMake](https://cmake.org/) - To build the project.
* A C++ compiler which supports C++20. Clang and GCC have be tested to work.



### Building

The CAMIRA code lives in a GitHub [repository](https://github.com/markageo/CAMIRA.git). To clone the repo, use the command:

    git clone https://github.com/markageo/CAMIRA.git

Then build CAMIRA like any other CMake project:

    mkdir build
    cd build
    cmake ../
    make 
Which will create an executable called `camira`. See [Usage](#usage) below for instructions on how to run the code. I recommend using `ccmake` or `cmake-gui` to configure so that the various compiler options can be seen and set interactively.


## Usage

After building the project, an executable called `camira` will be created. This can then be executed with an input file as:

    ./camira-flow flowInputFile.inp
    ./camira-plume plumeInputFile.inp

Where `flowInputFile.inp` and `plumeInputFile.inp` are the input files required by the flow and plume solvers respectively. All solver settings are controlled within these input files. 

A manual exists which describes the input file structure all the available options. This can be found in the `manual` directory of the [repository](https://github.com/markageo/camira.git). 

## Publications

Below is a list of publications related to the development of this project.

*  [George, M. A., Williamson, N., & Armfield, S. W. (2024). A coupled block implicit solver for the incompressible Navier–Stokes equations on collocated grids. Computers & Fluids, 284, 106426.](https://www.sciencedirect.com/science/article/pii/S0045793024002573)


* [George, M. A., Williamson, N., & Armfield, S. W. (2025). Mass-conserving ghost cell immersed boundary method with multigrid for coupled Navier-Stokes solvers. Journal of Computational Physics, 114276.](https://www.sciencedirect.com/science/article/pii/S0021999125005595)

* [George, M. A., Williamson, N., & Armfield, S. W. (2026). Coupled FAS multigrid for the incompressible Navier–Stokes equations on collocated grids. Computers & Fluids, 314, 107106.](https://www.sciencedirect.com/science/article/pii/S0045793026001489)


* [George, M. A., (2026) A Rapid Steady Solver for the Navier-Stokes Equations, Ph.D. thesis, The University of Sydney](https://ses.library.usyd.edu.au/handle/2123/35334)


## Roadmap

CAMIRA is a work in progress. Here are some features we would like to see appear in future:
* Inclusion of wall models for turbulence models.
* Addition of more sophisticated turbulence models, beyond the zero-equation models present. 
* Support for GPU parallelism, in both the flow and plume solvers.


## Contact

For any questions, ideas, suggestions, or bugs, please do not hesitate to contact me:

Mark George - **markgeorge0311@hotmail.com**



## Acknowledgements

My PhD Supervisors [Steven Armfield](https://www.sydney.edu.au/engineering/about/our-people/academic-staff/steven-armfield.html) and [Nicholas Williamson](https://www.sydney.edu.au/engineering/about/our-people/academic-staff/nicholas-williamson.html). Also [DMTC](https://dmtc.com.au/) and [DSTG Australia](https://www.dst.defence.gov.au/) for supporting the project.


