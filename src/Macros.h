#ifndef CFD_MACROS
#define CFD_MACROS


// Profiling macros
#ifdef CFD_PROFILING
#   if !defined(TIC) || !defined(TOC)
#       include "profiler/profiler.h"
#       define TIC(name) PROF::prof.tic(name);
#       define TOC(name) PROF::prof.toc(name);
#   endif

    namespace PROF {
        inline profiler<perf_counter::clock<time_units::SECONDS>> prof("CFD Profiling");
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
#   define CFD_PRAGMA_VECTORIZE _Pragma("clang loop vectorize(enable)")

# elif defined(__GNUC__) || defined(__GNUG__)
#   define CFD_PRAGMA_VECTORIZE _Pragma("GCC ivdep")

# else
#   define CFD_PRAGMA_VECTORIZE

# endif


// If compiled with -ffast-math (specifically -ffinite-math-only), compiler assumes nans cannot occur.
// This may be important in certain functions.
# if defined(__clang__)
#   if !defined(__FAST_MATH__)
#       define CFD_HONOR_INFINITIES_AND_NANS
#   endif

# elif defined(__GNUC__) || defined(__GNUG__)
#   if (HONOR_INFINITIES) && (HONOR_NANS)
#       define CFD_HONOR_INFINITIES_AND_NANS
#   endif

# else
#   define CFD_HONOR_INFINITIES_AND_NANS

# endif


#endif // CFD_MACROS