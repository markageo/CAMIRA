#ifndef PROFILER_H
#define PROFILER_H

#include <map>
#include <vector>
#include <string>
#include <optional>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "clock.h"

namespace PROF {
/// Profiler Class 
/**
 *  \param Counter  
 *  \param IDENT_level 
*/
template <class Counter = PROF::perf_counter::clock<> , size_t INDENT_LEVEL = 4>
class profiler {

public:
    typedef typename Counter::value_type value_type;
    typedef double                       delta_type;

    /// Initialize 
    profiler() : name("Profiler") {
        init();
    }

    profiler(const std::string& name) : name(name) {
        init();
    }

    void tic(const std::string& name){ 
        stack.back()->children[name].begin = counter.now();
        stack.back()->children[name].nCalls++;
        stack.push_back(&stack.back()->children[name]);
    }

    delta_type toc(const std::string& /*name*/ = ""){
        if (stack.size() == 1) return 0;

        auto top = stack.back();
        stack.pop_back();

        value_type now   = counter.now();
        delta_type delta = now - top->begin;

        top->time_delta += delta;
        root.time_delta  = now - root.begin;

        return delta;
    }

    void reset() {
        stack.clear();
        root.time_delta = 0;
        root.children.clear();

        stack.push_back(&root);
        root.begin = counter.now();
    }

private:
    struct profiler_unit {
        // Data 
        
        value_type begin;
        delta_type time_delta;
        size_t nCalls;
        std::map<std::string, profiler_unit> children;

        // functions 
        profiler_unit() : time_delta(0), nCalls(0) {}

        delta_type children_time() const {
            delta_type s = 0;
            for (const auto& [name, child]: children) 
                s += child.time_delta;
            return s;
        }

        void print(std::ostream& out, const std::string& name, int level, delta_type total, size_t width) const {
            using namespace fmt::literals;
            auto percentage_time = 100 * time_delta / total;
            constexpr auto format_string = "[{pad:>{level}}{name}:{pad:>{width}}{time:>15.{digits}f}{units}] [nCalls {counter}] ({percnt:>6.2f}%)\n"; 
            fmt::print(out, format_string,
                    "digits"_a = Counter::digits(),
                    "level"_a = level, 
                    "width"_a = width - level - name.size(), 
                    "units"_a = Counter::units(), 
                    "name"_a = name, 
                    "pad"_a = "", "time"_a = time_delta, 
                    "percnt"_a = percentage_time, "counter"_a = nCalls);

            if (children.size()) {
                delta_type val = time_delta - children_time();
                percentage_time = 100.0 * val / total;
                std::string str = "self";
                fmt::print(out, format_string,
                    "digits"_a = Counter::digits(),
                    "level"_a = level + INDENT_LEVEL/2, 
                    "width"_a = width - level - str.size() - INDENT_LEVEL/2, 
                    "units"_a = Counter::units(), 
                    "name"_a = str, 
                    "pad"_a = "", "time"_a = time_delta, 
                    "percnt"_a = percentage_time, "counter"_a = nCalls);
            }

            for (auto& [name, child] : children)
                child.print(out, name, level + INDENT_LEVEL, total, width);
        }

        size_t total_width(const std::string &name, int level) const {
            size_t w = name.size() + level;
            for(auto const& [name, child] : children)
                w = std::max(w, child.total_width(name, level + INDENT_LEVEL));
            return w;
        }
    };

    std::string name;
    Counter counter;
    profiler_unit root;
    std::vector<profiler_unit*> stack;

    void print(std::ostream &out) const {
        if (stack.back() != &root)
            fmt::print(out, "Warning! Profile is incomplete.\n");
        root.print(out, name, 0, root.time_delta, root.total_width(name, 0));
    }

    friend std::ostream& operator<<(std::ostream &out, const profiler &prof) {
        fmt::print(out, "\n");
        prof.print(out);
        return out;
    }
    
    void init() {
        stack.reserve(64);
        stack.push_back(&root);
        root.nCalls++;
        root.begin = counter.now();
    }
};

} // end namespace PROF

#endif