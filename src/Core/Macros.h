#ifndef CAMIRA_MACROS
#define CAMIRA_MACROS


// Profiling macros
#ifdef CAMIRA_PROFILING
#   if !defined(TIC) || !defined(TOC)
#       include "profiler/profiler.h"
#       define TIC(name) PROF::prof.tic(name);
#       define TOC(name) PROF::prof.toc(name);
#   endif

    namespace PROF {
        inline profiler<perf_counter::clock<time_units::SECONDS>> prof("CAMIRA_ Profiling");
    }

#else
#   ifndef TIC
#       define TIC(name)
#   endif
#   ifndef TOC
#       define TOC(name)
#   endif
#endif


// Loop autovectorisation
#if defined(__clang__)
#   define CAMIRA_PRAGMA_VECTORIZE _Pragma("clang loop vectorize(enable)")

# elif defined(__GNUC__) || defined(__GNUG__)
#   define CAMIRA_PRAGMA_VECTORIZE _Pragma("GCC ivdep")

# else
#   define CAMIRA_PRAGMA_VECTORIZE

# endif


// If compiled with -ffast-math (specifically -ffinite-math-only), compiler assumes nans cannot occur.
// This may be important in certain functions.
# if defined(__GNUC__) || defined(__GNUG__) || defined(__clang__)
#   if !defined(__FAST_MATH__)
#       define CAMIRA_HONOR_INFINITIES_AND_NANS
#   endif

# else
#   define CAMIRA_HONOR_INFINITIES_AND_NANS

# endif


#endif // CAMIRA_MACROS