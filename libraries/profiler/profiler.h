#ifndef PROFILER_H
#define PROFILER_H

#include <map>
#include <vector>
#include <string>
#include <optional>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "clock.h"

namespace PROF {
/// Profiler Class 
/**
 *  \param Counter  
*/
template <class Counter = PROF::perf_counter::clock<> >
class profiler {

public:
    typedef typename Counter::value_type value_type;
    typedef double                       delta_type;

    /// Initialize 
    profiler(const std::string& name, size_t indentLevel) : m_name(name), m_indentLevel(indentLevel) {
        init(); 
    }

    profiler() : profiler("Profiler", 4) {};

    profiler(const std::string& name) : profiler(name, 4) {};


    void tic(const std::string& name){ 
        m_stack.back()->children[name].begin = m_counter.now();
        m_stack.back()->children[name].nCalls++;
        m_stack.push_back(&m_stack.back()->children[name]);
    }

    delta_type toc(const std::string& /*name*/ = ""){
        if (m_stack.size() == 1) return 0;

        auto top = m_stack.back();
        m_stack.pop_back();

        value_type now   = m_counter.now();
        delta_type delta = now - top->begin;

        top->time_delta += delta;
        m_root.time_delta  = now - m_root.begin;

        return delta;
    }


    void reset() {
        m_stack.clear();
        m_root.time_delta = 0;
        m_root.children.clear();

        m_stack.push_back(&m_root);
        m_root.begin = m_counter.now();
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
            for (const auto& [childName, child]: children) 
                s += child.time_delta;
            return s;
        }

        void print(std::ostream& out, const std::string& name, int level, delta_type total, size_t nameWidth, size_t nCallsWidth, size_t indentLevel) const {
            
            auto formatted_line = [&]( std::string nameString, delta_type time_value, size_t indenting ) -> std::string {
                delta_type percentage_time = 100 * time_value / total;

                std::ostringstream formatted;
                formatted << "[" 
                          << std::setw(level + indenting) << std::right << ""                                         // Indenting
                          << nameString << ":"                                                                        // Name text
                          << std::setw(nameWidth - level - nameString.size() - indenting) << std::right << ""         // Padding
                          << std::setw(15) << std::fixed << std::setprecision(2) << time_value << Counter::units()    // Time
                          << "]"

                          << " [nCalls " << std::setw(nCallsWidth) << std::right  << nCalls << "]"               

                          << " ("        << std::setw(6) << std::right << std::fixed << std::setprecision(2) << percentage_time << "%)\n";

                return formatted.str();
            };


            // Total time
            out << formatted_line( name, time_delta, 0 );

            // Untimed sections
            if (!children.empty()) {
                delta_type delta_untimed = time_delta - children_time();
                out << formatted_line( "Untimed", delta_untimed, indentLevel );
            }

            // Timed sections
            for (const auto& [childName, child] : children)
                child.print(out, childName, level + indentLevel, total, nameWidth, nCallsWidth, indentLevel);
        }


        size_t total_width(const std::string &name, int level, size_t indentLevel) const {
            size_t w = name.size() + level;
            for(auto const& [childName, child] : children)
                w = std::max(w, child.total_width(childName, level + indentLevel, indentLevel));
            return w;
        }

        size_t num_digits( size_t num ) const {
            if ( num/10 == 0 )
                return 1;
            return 1 + num_digits(num / 10);
        }

        size_t nCalls_width() const {
            size_t w = num_digits( nCalls );
            for(auto const& [childName, child] : children)
                w = std::max( w, child.nCalls_width() );
            return w;
        }
    };

    std::string m_name;
    Counter m_counter;
    profiler_unit m_root;
    size_t m_indentLevel;
    std::vector<profiler_unit*> m_stack;

    void print(std::ostream &out) const {
        if (m_stack.back() != &m_root)
            out << "Warning! Profile is incomplete.\n";
        m_root.print(out, m_name, 0, m_root.time_delta, m_root.total_width(m_name, 0, m_indentLevel), m_root.nCalls_width(), m_indentLevel);
    }

    friend std::ostream& operator<<(std::ostream &out, const profiler &prof) {
        out << "\n";
        prof.print(out);
        return out;
    }
    
    void init() {
        m_stack.reserve(64);
        m_stack.push_back(&m_root);
        m_root.nCalls++;
        m_root.begin = m_counter.now();
    }
};

} // end namespace PROF

#endif