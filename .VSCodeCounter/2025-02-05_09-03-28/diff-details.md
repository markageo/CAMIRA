# Diff Details

Date : 2025-02-05 09:03:28

Directory /home/mark/PhD/Code/coupled-sweep-ns/src

Total : 113 files,  239 codes, 17 comments, 87 blanks, all 343 lines

[Summary](results.md) / [Details](details.md) / [Diff Summary](diff.md) / Diff Details

## Files
| filename | language | code | comment | blank | total |
| :--- | :--- | ---: | ---: | ---: | ---: |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/CMakeLists.txt](//home/mark/PhD/Code/coupled-sweep-ns-main/src/CMakeLists.txt) | CMake | -23 | 0 | -9 | -32 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/CMakeLists.txt](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/CMakeLists.txt) | CMake | -13 | 0 | 0 | -13 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FaceInterpolatedVelocity.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FaceInterpolatedVelocity.h) | C++ | -89 | -4 | -33 | -126 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FaceVelocities.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FaceVelocities.cpp) | C++ | -301 | -25 | -110 | -436 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolume.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolume.h) | C++ | -5 | 0 | -3 | -8 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolumeCoefficients.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolumeCoefficients.cpp) | C++ | -1,681 | -209 | -684 | -2,574 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolumeFunctions.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolumeFunctions.h) | C++ | -68 | -11 | -32 | -111 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolumeStructures.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolumeStructures.cpp) | C++ | -341 | -24 | -118 | -483 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolumeStructures.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/FiniteVolumeStructures.h) | C++ | -61 | -5 | -15 | -81 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/GhostCells.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/GhostCells.cpp) | C++ | -60 | -2 | -22 | -84 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/Mesh.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/Mesh.cpp) | C++ | -376 | -23 | -183 | -582 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/Mesh.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/Mesh.h) | C++ | -33 | -4 | -11 | -48 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/VertexValues.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/FiniteVolume/VertexValues.cpp) | C++ | -208 | -22 | -88 | -318 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Geometry/CMakeLists.txt](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Geometry/CMakeLists.txt) | CMake | -5 | 0 | 0 | -5 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Geometry/Geometry.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Geometry/Geometry.cpp) | C++ | -223 | -42 | -106 | -371 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Geometry/Geometry.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Geometry/Geometry.h) | C++ | -28 | 0 | -17 | -45 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/ArrayIO.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/ArrayIO.h) | C++ | -101 | -8 | -39 | -148 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/CMakeLists.txt](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/CMakeLists.txt) | CMake | -13 | 0 | 0 | -13 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/CSVReader.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/CSVReader.cpp) | C++ | -96 | 0 | -33 | -129 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/CSVReader.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/CSVReader.h) | C++ | -6 | 0 | -3 | -9 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/IOTools.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/IOTools.h) | C++ | -41 | -6 | -19 | -66 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/InputParser.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/InputParser.cpp) | C++ | -348 | -49 | -117 | -514 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/InputParser.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/InputParser.h) | C++ | -10 | -1 | -7 | -18 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/InputProcessing.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/InputProcessing.cpp) | C++ | -665 | -94 | -236 | -995 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/InputProcessing.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/InputProcessing.h) | C++ | -104 | -11 | -36 | -151 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/VTKReader.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/VTKReader.cpp) | C++ | -113 | -8 | -59 | -180 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/VTKReader.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/VTKReader.h) | C++ | -15 | 0 | -12 | -27 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/VTKWriter.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/IO/VTKWriter.h) | C++ | -286 | -33 | -96 | -415 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/ImmersedBoundary/CMakeLists.txt](//home/mark/PhD/Code/coupled-sweep-ns-main/src/ImmersedBoundary/CMakeLists.txt) | CMake | -6 | 0 | 0 | -6 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/ImmersedBoundary/IBData.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/ImmersedBoundary/IBData.cpp) | C++ | -237 | -20 | -109 | -366 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/ImmersedBoundary/IBSolverFunctions.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/ImmersedBoundary/IBSolverFunctions.cpp) | C++ | -101 | -8 | -57 | -166 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/ImmersedBoundary/ImmersedBoundary.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/ImmersedBoundary/ImmersedBoundary.h) | C++ | -49 | -12 | -34 | -95 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Macros.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Macros.h) | C++ | -34 | -4 | -15 | -53 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Multigrid/CMakeLists.txt](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Multigrid/CMakeLists.txt) | CMake | -6 | 0 | 0 | -6 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Multigrid/Multigrid.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Multigrid/Multigrid.cpp) | C++ | -457 | -58 | -175 | -690 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Multigrid/Multigrid.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Multigrid/Multigrid.h) | C++ | -65 | -1 | -22 | -88 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/CMakeLists.txt](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/CMakeLists.txt) | CMake | -13 | 0 | 0 | -13 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/ConvergenceLogging.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/ConvergenceLogging.h) | C++ | -302 | -13 | -88 | -403 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/LineSolver.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/LineSolver.h) | C++ | -95 | -6 | -27 | -128 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/LinearSolver.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/LinearSolver.h) | C++ | -227 | -17 | -81 | -325 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/PlaneSolver.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/PlaneSolver.h) | C++ | -83 | -4 | -30 | -117 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/ResidualFunctions.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/ResidualFunctions.h) | C++ | -285 | -25 | -108 | -418 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/Solver.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/Solver.cpp) | C++ | -455 | -33 | -124 | -612 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/Solver.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/Solver.h) | C++ | -14 | 0 | -7 | -21 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/StaggerIndexing.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/StaggerIndexing.h) | C++ | -117 | -2 | -25 | -144 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/StencilConstants.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/StencilConstants.h) | C++ | -158 | -14 | -59 | -231 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/TriadSolver.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Solver/TriadSolver.h) | C++ | -283 | -29 | -106 | -418 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/CMakeLists.txt](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/CMakeLists.txt) | CMake | -8 | 0 | 0 | -8 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/FVLookups.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/FVLookups.h) | C++ | -53 | -7 | -23 | -83 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/FVTools.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/FVTools.h) | C++ | -56 | -11 | -29 | -96 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/FieldProbe.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/FieldProbe.h) | C++ | -105 | -17 | -43 | -165 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/SweepTransformations.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/SweepTransformations.cpp) | C++ | -438 | -64 | -223 | -725 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/SweepTransformations.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Tools/SweepTransformations.h) | C++ | -38 | -18 | -29 | -85 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/Types.h](//home/mark/PhD/Code/coupled-sweep-ns-main/src/Types.h) | C++ | -237 | -27 | -88 | -352 |
| [/home/mark/PhD/Code/coupled-sweep-ns-main/src/main.cpp](//home/mark/PhD/Code/coupled-sweep-ns-main/src/main.cpp) | C++ | -73 | -15 | -24 | -112 |
| [src/CMakeLists.txt](/src/CMakeLists.txt) | CMake | 24 | 0 | 9 | 33 |
| [src/DerivedQuantities/CMakeLists.txt](/src/DerivedQuantities/CMakeLists.txt) | CMake | 6 | 0 | 0 | 6 |
| [src/DerivedQuantities/DerivedQuantities.h](/src/DerivedQuantities/DerivedQuantities.h) | C++ | 37 | 2 | 13 | 52 |
| [src/DerivedQuantities/FieldProbe.cpp](/src/DerivedQuantities/FieldProbe.cpp) | C++ | 94 | 16 | 42 | 152 |
| [src/DerivedQuantities/ForceCalculator.cpp](/src/DerivedQuantities/ForceCalculator.cpp) | C++ | 131 | 6 | 48 | 185 |
| [src/FiniteVolume/CMakeLists.txt](/src/FiniteVolume/CMakeLists.txt) | CMake | 13 | 0 | 0 | 13 |
| [src/FiniteVolume/FaceInterpolatedVelocity.h](/src/FiniteVolume/FaceInterpolatedVelocity.h) | C++ | 89 | 4 | 33 | 126 |
| [src/FiniteVolume/FaceVelocities.cpp](/src/FiniteVolume/FaceVelocities.cpp) | C++ | 301 | 25 | 110 | 436 |
| [src/FiniteVolume/FiniteVolume.h](/src/FiniteVolume/FiniteVolume.h) | C++ | 5 | 0 | 3 | 8 |
| [src/FiniteVolume/FiniteVolumeCoefficients.cpp](/src/FiniteVolume/FiniteVolumeCoefficients.cpp) | C++ | 1,678 | 209 | 684 | 2,571 |
| [src/FiniteVolume/FiniteVolumeFunctions.h](/src/FiniteVolume/FiniteVolumeFunctions.h) | C++ | 68 | 11 | 32 | 111 |
| [src/FiniteVolume/FiniteVolumeStructures.cpp](/src/FiniteVolume/FiniteVolumeStructures.cpp) | C++ | 341 | 24 | 118 | 483 |
| [src/FiniteVolume/FiniteVolumeStructures.h](/src/FiniteVolume/FiniteVolumeStructures.h) | C++ | 61 | 5 | 15 | 81 |
| [src/FiniteVolume/GhostCells.cpp](/src/FiniteVolume/GhostCells.cpp) | C++ | 60 | 2 | 22 | 84 |
| [src/FiniteVolume/Mesh.cpp](/src/FiniteVolume/Mesh.cpp) | C++ | 376 | 23 | 183 | 582 |
| [src/FiniteVolume/Mesh.h](/src/FiniteVolume/Mesh.h) | C++ | 33 | 4 | 11 | 48 |
| [src/FiniteVolume/VertexValues.cpp](/src/FiniteVolume/VertexValues.cpp) | C++ | 208 | 22 | 88 | 318 |
| [src/Geometry/CMakeLists.txt](/src/Geometry/CMakeLists.txt) | CMake | 5 | 0 | 0 | 5 |
| [src/Geometry/Geometry.cpp](/src/Geometry/Geometry.cpp) | C++ | 223 | 42 | 106 | 371 |
| [src/Geometry/Geometry.h](/src/Geometry/Geometry.h) | C++ | 28 | 0 | 17 | 45 |
| [src/IO/ArrayIO.h](/src/IO/ArrayIO.h) | C++ | 101 | 8 | 39 | 148 |
| [src/IO/CMakeLists.txt](/src/IO/CMakeLists.txt) | CMake | 13 | 0 | 0 | 13 |
| [src/IO/CSVReader.cpp](/src/IO/CSVReader.cpp) | C++ | 96 | 0 | 33 | 129 |
| [src/IO/CSVReader.h](/src/IO/CSVReader.h) | C++ | 6 | 0 | 3 | 9 |
| [src/IO/IOTools.h](/src/IO/IOTools.h) | C++ | 41 | 6 | 19 | 66 |
| [src/IO/InputParser.cpp](/src/IO/InputParser.cpp) | C++ | 348 | 49 | 117 | 514 |
| [src/IO/InputParser.h](/src/IO/InputParser.h) | C++ | 10 | 1 | 7 | 18 |
| [src/IO/InputProcessing.cpp](/src/IO/InputProcessing.cpp) | C++ | 676 | 95 | 241 | 1,012 |
| [src/IO/InputProcessing.h](/src/IO/InputProcessing.h) | C++ | 107 | 11 | 37 | 155 |
| [src/IO/VTKReader.cpp](/src/IO/VTKReader.cpp) | C++ | 113 | 8 | 59 | 180 |
| [src/IO/VTKReader.h](/src/IO/VTKReader.h) | C++ | 15 | 0 | 12 | 27 |
| [src/IO/VTKWriter.h](/src/IO/VTKWriter.h) | C++ | 286 | 33 | 96 | 415 |
| [src/ImmersedBoundary/CMakeLists.txt](/src/ImmersedBoundary/CMakeLists.txt) | CMake | 6 | 0 | 0 | 6 |
| [src/ImmersedBoundary/IBData.cpp](/src/ImmersedBoundary/IBData.cpp) | C++ | 252 | 22 | 111 | 385 |
| [src/ImmersedBoundary/IBSolverFunctions.cpp](/src/ImmersedBoundary/IBSolverFunctions.cpp) | C++ | 101 | 8 | 57 | 166 |
| [src/ImmersedBoundary/ImmersedBoundary.h](/src/ImmersedBoundary/ImmersedBoundary.h) | C++ | 49 | 12 | 34 | 95 |
| [src/Macros.h](/src/Macros.h) | C++ | 34 | 4 | 15 | 53 |
| [src/Multigrid/CMakeLists.txt](/src/Multigrid/CMakeLists.txt) | CMake | 6 | 0 | 0 | 6 |
| [src/Multigrid/Multigrid.cpp](/src/Multigrid/Multigrid.cpp) | C++ | 457 | 58 | 175 | 690 |
| [src/Multigrid/Multigrid.h](/src/Multigrid/Multigrid.h) | C++ | 65 | 1 | 22 | 88 |
| [src/Solver/CMakeLists.txt](/src/Solver/CMakeLists.txt) | CMake | 13 | 0 | 0 | 13 |
| [src/Solver/ConvergenceLogging.h](/src/Solver/ConvergenceLogging.h) | C++ | 332 | 15 | 98 | 445 |
| [src/Solver/LineSolver.h](/src/Solver/LineSolver.h) | C++ | 95 | 6 | 27 | 128 |
| [src/Solver/LinearSolver.h](/src/Solver/LinearSolver.h) | C++ | 227 | 17 | 84 | 328 |
| [src/Solver/PlaneSolver.h](/src/Solver/PlaneSolver.h) | C++ | 83 | 4 | 30 | 117 |
| [src/Solver/ResidualFunctions.h](/src/Solver/ResidualFunctions.h) | C++ | 285 | 25 | 108 | 418 |
| [src/Solver/Solver.cpp](/src/Solver/Solver.cpp) | C++ | 472 | 38 | 129 | 639 |
| [src/Solver/Solver.h](/src/Solver/Solver.h) | C++ | 14 | 0 | 7 | 21 |
| [src/Solver/StaggerIndexing.h](/src/Solver/StaggerIndexing.h) | C++ | 117 | 2 | 25 | 144 |
| [src/Solver/StencilConstants.h](/src/Solver/StencilConstants.h) | C++ | 158 | 14 | 59 | 231 |
| [src/Solver/TriadSolver.h](/src/Solver/TriadSolver.h) | C++ | 283 | 29 | 106 | 418 |
| [src/Tools/CMakeLists.txt](/src/Tools/CMakeLists.txt) | CMake | 7 | 0 | 0 | 7 |
| [src/Tools/FVLookups.h](/src/Tools/FVLookups.h) | C++ | 53 | 7 | 23 | 83 |
| [src/Tools/FVTools.h](/src/Tools/FVTools.h) | C++ | 56 | 11 | 29 | 96 |
| [src/Tools/SweepTransformations.cpp](/src/Tools/SweepTransformations.cpp) | C++ | 438 | 64 | 223 | 725 |
| [src/Tools/SweepTransformations.h](/src/Tools/SweepTransformations.h) | C++ | 38 | 18 | 29 | 85 |
| [src/Types.h](/src/Types.h) | C++ | 240 | 27 | 89 | 356 |
| [src/main.cpp](/src/main.cpp) | C++ | 73 | 15 | 24 | 112 |

[Summary](results.md) / [Details](details.md) / [Diff Summary](diff.md) / Diff Details