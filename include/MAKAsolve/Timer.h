#ifndef MAKASOLVE_TIMER_H
#define MAKASOLVE_TIMER_H

#include <PCU.h>
#include <chrono>
#include <list>
#include <sstream>
#include <string>

typedef unsigned long long ticks;

namespace maka {

class Timer {
public:
	Timer(pcu::PCU* pcu, int precision = 9);

	void start_time();

	void stop_time(std::string header);

	// prints times + prefix data in a single line
	void print_times_line();

	// prints headers
	void print_header_line();

	// prepend some info (e.g. to associate times with a run)
	template <typename T> void prepend_info(std::string header, T data) {
		header_prefix_ += header + delim_;

		std::ostringstream oss;
		oss << prefix_ << data << delim_;
		prefix_ = oss.str();
	}

private:
	pcu::PCU* pcu_;

	std::list<std::string> headers_;
	std::list<double> time_ms_;
	std::string delim_ = " ";

	std::string header_prefix_;
	std::string prefix_;

	std::chrono::steady_clock::time_point last_time_;
};

} // namespace maka

#endif // MAKASOLVE_TIMER_H