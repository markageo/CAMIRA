#ifndef CFD_MACROS
#define CFD_MACROS


// Profiling macros
#ifdef PROFILING
#   if !defined(TIC) || !defined(TOC)
#       include "profiler/profiler.h"
#       define TIC(name) PROF::prof.tic(name);
#       define TOC(name) PROF::prof.toc(name);
#   endif

    namespace PROF {
        inline profiler<perf_counter::clock<time_units::SECONDS>> prof;
    }

#else
#   ifndef TIC
#       define TIC(name)
#   endif
#   ifndef TOC
#       define TOC(name)
#   endif
#endif


// Compiler specific macros
#if defined(__clang__)
#   define CFD_PRAGMA_VECTORIZE _Pragma("clang loop vectorize(enable)")
# elif defined(__GNUC__) || defined(__GNUG__)
#   define CFD_PRAGMA_VECTORIZE _Pragma("GCC ivdep")
# else
#   define CFD_PRAGMA_VECTORIZE
# endif



#endif // CFD_MACROS