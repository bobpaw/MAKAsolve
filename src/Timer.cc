#include <MAKAsolve/Timer.h>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace maka {

Timer::Timer(pcu::PCU* pcu, int precision) : pcu_(pcu) {
	std::cout << std::setprecision(precision);
}

void Timer::start_time() {
	pcu_->Barrier();
	last_time_ = std::chrono::steady_clock::now();
}

void Timer::stop_time(std::string header) {
	pcu_->Barrier();
	auto now_time = std::chrono::steady_clock::now();
	std::chrono::duration<double, std::milli> elapsed = now_time - last_time_;

	headers_.push_back(header);
	time_ms_.push_back(elapsed.count());
}

// prints times + prefix data in a single line
void Timer::print_times_line() {
	std::cout << prefix_;
	for (double t : time_ms_) {
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