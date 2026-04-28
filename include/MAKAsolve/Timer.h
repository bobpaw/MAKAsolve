#ifndef MAKASOLVE_TIMER_H
#define MAKASOLVE_TIMER_H

#include <list>
#include <string>
#include <PCU.h>
#include <sstream>

typedef unsigned long long ticks;

namespace maka {

class Timer {
public:
    Timer(int precision = 9);

    void start_time();

    void stop_time(std::string header);

    // prints times in a single line
    void print_times_line();

    // prints headers
    void print_header_line();

    // prepend some info (e.g. to associate times with a run)
    template <typename T>
    void prepend_info(std::string header, T data) {
        header_prefix_ += header + delim_;

        std::ostringstream oss;
        oss << prefix_ << data << delim_;
        prefix_ = oss.str();
    }

private:
    std::list<std::string> headers_;
    std::list<double> time_s_;
    std::string delim_ = " ";

    std::string header_prefix_;
    std::string prefix_;

    ticks start_ticks_;
    
};

} // namespace maka

#endif // MAKASOLVE_TIMER_H