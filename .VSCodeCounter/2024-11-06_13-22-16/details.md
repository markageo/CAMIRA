# Details

Date : 2024-11-06 13:22:16

Directory /home/mark/PhD/Code/coupled-sweep-ns-patterns/src

Total : 55 files,  9070 codes, 986 comments, 3545 blanks, all 13601 lines

[Summary](results.md) / Details / [Diff Summary](diff.md) / [Diff Details](diff-details.md)

## Files
| filename | language | code | comment | blank | total |
| :--- | :--- | ---: | ---: | ---: | ---: |
| [src/CMakeLists.txt](/src/CMakeLists.txt) | CMake | 23 | 0 | 9 | 32 |
| [src/FiniteVolume/CMakeLists.txt](/src/FiniteVolume/CMakeLists.txt) | CMake | 13 | 0 | 0 | 13 |
| [src/FiniteVolume/FaceInterpolatedVelocity.h](/src/FiniteVolume/FaceInterpolatedVelocity.h) | C++ | 89 | 4 | 33 | 126 |
| [src/FiniteVolume/FaceVelocities.cpp](/src/FiniteVolume/FaceVelocities.cpp) | C++ | 290 | 25 | 108 | 423 |
| [src/FiniteVolume/FiniteVolume.h](/src/FiniteVolume/FiniteVolume.h) | C++ | 5 | 0 | 3 | 8 |
| [src/FiniteVolume/FiniteVolumeCoefficients.cpp](/src/FiniteVolume/FiniteVolumeCoefficients.cpp) | C++ | 1,583 | 205 | 665 | 2,453 |
| [src/FiniteVolume/FiniteVolumeFunctions.h](/src/FiniteVolume/FiniteVolumeFunctions.h) | C++ | 68 | 11 | 32 | 111 |
| [src/FiniteVolume/FiniteVolumeStructures.cpp](/src/FiniteVolume/FiniteVolumeStructures.cpp) | C++ | 341 | 24 | 118 | 483 |
| [src/FiniteVolume/FiniteVolumeStructures.h](/src/FiniteVolume/FiniteVolumeStructures.h) | C++ | 60 | 5 | 15 | 80 |
| [src/FiniteVolume/GhostCells.cpp](/src/FiniteVolume/GhostCells.cpp) | C++ | 46 | 1 | 21 | 68 |
| [src/FiniteVolume/Mesh.cpp](/src/FiniteVolume/Mesh.cpp) | C++ | 371 | 21 | 181 | 573 |
| [src/FiniteVolume/Mesh.h](/src/FiniteVolume/Mesh.h) | C++ | 33 | 4 | 11 | 48 |
| [src/FiniteVolume/VertexValues.cpp](/src/FiniteVolume/VertexValues.cpp) | C++ | 197 | 22 | 84 | 303 |
| [src/Geometry/CMakeLists.txt](/src/Geometry/CMakeLists.txt) | CMake | 5 | 0 | 0 | 5 |
| [src/Geometry/Geometry.cpp](/src/Geometry/Geometry.cpp) | C++ | 163 | 30 | 76 | 269 |
| [src/Geometry/Geometry.h](/src/Geometry/Geometry.h) | C++ | 13 | 1 | 9 | 23 |
| [src/IO/ArrayIO.h](/src/IO/ArrayIO.h) | C++ | 101 | 8 | 39 | 148 |
| [src/IO/CMakeLists.txt](/src/IO/CMakeLists.txt) | CMake | 13 | 0 | 0 | 13 |
| [src/IO/CSVReader.cpp](/src/IO/CSVReader.cpp) | C++ | 96 | 0 | 33 | 129 |
| [src/IO/CSVReader.h](/src/IO/CSVReader.h) | C++ | 6 | 0 | 3 | 9 |
| [src/IO/IOTools.h](/src/IO/IOTools.h) | C++ | 41 | 6 | 19 | 66 |
| [src/IO/InputParser.cpp](/src/IO/InputParser.cpp) | C++ | 348 | 49 | 117 | 514 |
| [src/IO/InputParser.h](/src/IO/InputParser.h) | C++ | 10 | 1 | 7 | 18 |
| [src/IO/InputProcessing.cpp](/src/IO/InputProcessing.cpp) | C++ | 629 | 93 | 221 | 943 |
| [src/IO/InputProcessing.h](/src/IO/InputProcessing.h) | C++ | 105 | 11 | 34 | 150 |
| [src/IO/VTKReader.cpp](/src/IO/VTKReader.cpp) | C++ | 113 | 8 | 59 | 180 |
| [src/IO/VTKReader.h](/src/IO/VTKReader.h) | C++ | 15 | 0 | 12 | 27 |
| [src/IO/VTKWriter.h](/src/IO/VTKWriter.h) | C++ | 286 | 33 | 96 | 415 |
| [src/ImmersedBoundary/CMakeLists.txt](/src/ImmersedBoundary/CMakeLists.txt) | CMake | 6 | 0 | 0 | 6 |
| [src/ImmersedBoundary/IBData.cpp](/src/ImmersedBoundary/IBData.cpp) | C++ | 265 | 24 | 119 | 408 |
| [src/ImmersedBoundary/IBSolverFunctions.cpp](/src/ImmersedBoundary/IBSolverFunctions.cpp) | C++ | 101 | 8 | 57 | 166 |
| [src/ImmersedBoundary/ImmersedBoundary.h](/src/ImmersedBoundary/ImmersedBoundary.h) | C++ | 48 | 12 | 34 | 94 |
| [src/Macros.h](/src/Macros.h) | C++ | 38 | 4 | 16 | 58 |
| [src/Multigrid/CMakeLists.txt](/src/Multigrid/CMakeLists.txt) | CMake | 6 | 0 | 0 | 6 |
| [src/Multigrid/Multigrid.cpp](/src/Multigrid/Multigrid.cpp) | C++ | 450 | 58 | 173 | 681 |
| [src/Multigrid/Multigrid.h](/src/Multigrid/Multigrid.h) | C++ | 65 | 1 | 22 | 88 |
| [src/Solver/CMakeLists.txt](/src/Solver/CMakeLists.txt) | CMake | 13 | 0 | 0 | 13 |
| [src/Solver/ConvergenceLogging.h](/src/Solver/ConvergenceLogging.h) | C++ | 293 | 10 | 86 | 389 |
| [src/Solver/LineSolver.h](/src/Solver/LineSolver.h) | C++ | 95 | 6 | 27 | 128 |
| [src/Solver/LinearSolver.h](/src/Solver/LinearSolver.h) | C++ | 241 | 20 | 87 | 348 |
| [src/Solver/PlaneSolver.h](/src/Solver/PlaneSolver.h) | C++ | 83 | 4 | 30 | 117 |
| [src/Solver/ResidualFunctions.h](/src/Solver/ResidualFunctions.h) | C++ | 292 | 26 | 109 | 427 |
| [src/Solver/Solver.cpp](/src/Solver/Solver.cpp) | C++ | 458 | 37 | 125 | 620 |
| [src/Solver/Solver.h](/src/Solver/Solver.h) | C++ | 14 | 0 | 7 | 21 |
| [src/Solver/StaggerIndexing.h](/src/Solver/StaggerIndexing.h) | C++ | 117 | 2 | 25 | 144 |
| [src/Solver/StencilConstants.h](/src/Solver/StencilConstants.h) | C++ | 158 | 14 | 59 | 231 |
| [src/Solver/TriadSolver.h](/src/Solver/TriadSolver.h) | C++ | 272 | 39 | 111 | 422 |
| [src/Tools/CMakeLists.txt](/src/Tools/CMakeLists.txt) | CMake | 8 | 0 | 0 | 8 |
| [src/Tools/FVLookups.h](/src/Tools/FVLookups.h) | C++ | 53 | 7 | 23 | 83 |
| [src/Tools/FVTools.h](/src/Tools/FVTools.h) | C++ | 56 | 11 | 29 | 96 |
| [src/Tools/FieldProbe.h](/src/Tools/FieldProbe.h) | C++ | 105 | 17 | 43 | 165 |
| [src/Tools/SweepTransformations.cpp](/src/Tools/SweepTransformations.cpp) | C++ | 427 | 65 | 221 | 713 |
| [src/Tools/SweepTransformations.h](/src/Tools/SweepTransformations.h) | C++ | 36 | 17 | 26 | 79 |
| [src/Types.h](/src/Types.h) | C++ | 236 | 27 | 88 | 351 |
| [src/main.cpp](/src/main.cpp) | C++ | 71 | 15 | 23 | 109 |

[Summary](results.md) / Details / [Diff Summary](diff.md) / [Diff Details](diff-details.md)