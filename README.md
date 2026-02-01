
<!-- # <img src="images/logo.png?raw=true" width="128" valign="middle" alt="RAJA"/> -->

# CAMIRA
CAMIRA (**C**oupled **A**lgorithm for **M**ultigrid **I**mmersed Boundary **R**eynolds **A**veraged Navier-Stokes, also an Aboriginal name for the word "windy") is a highly efficiency solver for the incompressible steady-state Reynolds-Averaged Navier-Stokes equations. The code and the methods used were developed as part of my PhD. 

CAMIRA uses a fully coupled matrix free smoother with Full Approximation Scheme (FAS) multigrid to solve the equations on a rectillinear collocated finite volume grid. Complex geometries are accounted for using a mass conservative directional immersed boundary method. Details on these numerical methods can be found in the various [publications](#publications) that are part of this work. CAMIRA can run in parallel on shared memory CPUs using OpenMP.


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

    ./camira inputFile.inp

Where `inputFile.inp` is the input file used to run the code. All solver settings are controlled through this single input file. 

A manual exists which describes how to create an input file and all the available options. This can be found in the `manual` directory of the [repository](https://github.com/mgeo2280/camira.git). 

## Publications

Below is a list of publications related to the development of this project.

*  [George, M. A., Williamson, N., & Armfield, S. W. (2024). A coupled block implicit solver for the incompressible Navier–Stokes equations on collocated grids. Computers & Fluids, 284, 106426.](https://www.sciencedirect.com/science/article/pii/S0045793024002573)


* [George, M. A., Williamson, N., & Armfield, S. W. (2025). Mass-conserving ghost cell immersed boundary method with multigrid for coupled Navier-Stokes solvers. Journal of Computational Physics, 114276.](https://www.sciencedirect.com/science/article/pii/S0021999125005595)

* [George, M., Williamson, N., & Armfield, S. W. Coupled FAS Multigrid for the Incompressible Navier-Stokes Equations on Collocated Grids. Available at SSRN 5460400 (Under review in Computers & Fluids).](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5460400)


* [George, M. A., (2026) A Rapid Steady Solver for the Navier-Stokes Equations, Ph.D. thesis, The University of Sydney]()


## Roadmap

CAMIRA is a work in progress. Here are some upcoming features:
* Inclusion of wall models for turbulence models.
* Addition of more sophisticated turbulence models, beyond the zero-equation models present. 
* Support for GPU and CPU parallelism through the [RAJA Portability Suite](https://github.com/LLNL/RAJA).


## Contact

For any questions, ideas, suggestions, or bugs, please do not hesitate to contact me:

Mark George - **markgeorge0311@hotmail.com**



## Acknowledgements

My PhD Supervisors [Steven Armfield](https://www.sydney.edu.au/engineering/about/our-people/academic-staff/steven-armfield.html) and [Nicholas Williamson](https://www.sydney.edu.au/engineering/about/our-people/academic-staff/nicholas-williamson.html). Also [DMTC](https://dmtc.com.au/) and [DSTG Australia](https://www.dst.defence.gov.au/) for supporting the project.


License
-----------

Distributed under the project_license. See `LICENSE.txt` for more information.