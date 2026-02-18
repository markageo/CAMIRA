target_sources( libs 
    PRIVATE 
        Solver.cpp 
    PRIVATE 
        Solver.h 
        LinearSolver.h 
        PlaneSolver.h 
        LineSolver.h 
        TriadSolver.h 
        StencilConstants.h
        ResidualFunctions.h 
        StaggerIndexing.h 
        ConvergenceLogging.h )