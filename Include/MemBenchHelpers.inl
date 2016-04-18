
template<const int M> int  mem_bench_impl::MemBenchTimers::check_tick() {
	
	
	int minDelta{};
	int Delta{};
	double t1{};
	double t2{};
	double timesfound[M] = { 0.0 };
	constexpr double EPS{ 1.0E-6 };
	/* Collects now a sequence of M unique values from the system timer */
	for (int i{ 0 }; i != M; ++i) {

		t1 = timer_tick();
		while ((std::abs(t2 = timer_tick() - t1)) <= EPS)
			;
		timesfound[i] = t1 = t2;
	}
	minDelta = 1000000;
	for (int i{ 1 }; i != M; ++i) {
		Delta = static_cast<int>(1.0E6 * (timesfound[i] - timesfound[i - 1]));
		minDelta = std::min<int>(minDelta, std::max<int>(Delta, 0));
	}
#ifndef MEM_BENCH_VERBOSE
	std::cout << "Execting in: MemBenchTimers::check_tick() " << std::endl;
#endif
	return (minDelta);
}

auto  mem_bench_impl::MemBenchTimers::timer_tick()->double {
	
	LARGE_INTEGER freq, time;
	::QueryPerformanceFrequency(&freq);
	::QueryPerformanceCounter(&time);
#ifndef MEM_BENCH_VERBOSE
	std::cout << "Executing in: MemBenchTimers::timer_tick() " << std::endl;
#endif
	return (time.QuadPart / static_cast<double>(freq.QuadPart));
}