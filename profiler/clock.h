#ifndef COUNTER_CLOCK_H
#define COUNTER_CLOCK_H

#include <ctime>
#include <chrono>

#include <fmt/core.h>

namespace PROF {

enum time_units {
    SECONDS,
    MILLISECONDS,
    NANOSECONDS
};

namespace perf_counter {

    template<time_units unit = time_units::SECONDS>
    struct clock {
        typedef double value_type;

        // static const char* units() { return "s"; }
        static const char* units() {
            if constexpr (unit == time_units::SECONDS)
                return "s";
            else if constexpr (unit == time_units::MILLISECONDS)
                return "ms";
            else if constexpr (unit == time_units::NANOSECONDS)
                return "ns";
        }

        static int digits() {
            if constexpr (unit == time_units::SECONDS)
                return 3;
            else if constexpr (unit == time_units::MILLISECONDS)
                return 0;
            else if constexpr (unit == time_units::NANOSECONDS)
                return 0;
        }

        static double now() {
            using namespace std::chrono;
            auto t = high_resolution_clock::now();
            
            if constexpr (unit == time_units::SECONDS) {
                auto delta = duration<double>(t.time_since_epoch());
                return delta.count();
            }
            else if constexpr (unit == time_units::MILLISECONDS) {
                auto delta = duration_cast<milliseconds>(t.time_since_epoch());
                return delta.count();
            }
            else if constexpr (unit == time_units::NANOSECONDS) {
                auto delta = duration_cast<nanoseconds>(t.time_since_epoch());
                return delta.count();
            }
        }
    };

} // end namespace perf_counter 
} // end namespace PROF


#endif // end COUNTER_CLOCK_H