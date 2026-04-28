#include <MAKAsolve/Timer.h>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace maka {

#ifdef POWER9
// IBM POWER9 System clock with 512MHZ resolution.
static __inline__ ticks getticks(void) {
	unsigned int tbl, tbu0, tbu1;

	do {
		__asm__ __volatile__("mftbu %0" : "=r"(tbu0));
		__asm__ __volatile__("mftb %0" : "=r"(tbl));
		__asm__ __volatile__("mftbu %0" : "=r"(tbu1));
	} while (tbu0 != tbu1);

	return (((unsigned long long)tbu0) << 32) | tbl;
}
#endif

Timer::Timer(int precision) {
	std::cout << std::setprecision(precision);
}

void Timer::start_time() {
#ifdef POWER9
	start_ticks_ = getticks();
#else
// if we want to time outside of AiMOS
#endif
}

void Timer::stop_time(std::string header) {
	double time = 0;
#ifdef POWER9
	ticks stop_ticks_ = getticks();
	time = (double)(stop_ticks_ - start_ticks_) / (double)512000000.0;
#else
// if we want to time outside of AiMOS
#endif

	headers_.push_back(header);
	time_s_.push_back(time);
}

// prints times + prefix data in a single line
void Timer::print_times_line() {
	std::cout << prefix_;
	for (double t : time_s_) {
		std::cout << t << delim_;
	}
	std::cout << std::endl;
}

// prints headers
void Timer::print_header_line() {
	std::cout << header_prefix_;
	for (std::string s : headers_) {
		std::cout << s << delim_;
	}
	std::cout << std::endl;
}

} // namespace maka