
#include "MemBenchAlignMalloc.h"

/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Memory Benchmarking tool adapted from John D. McCalpin "STREAM" benchmark
@file MemBenchNewPtr.cpp
@author: Bernard Gingold
@version:  1.0  15/04/2016 19:59
@description: MemBenchAlignMalloc.h
@reference: STREAM Benchmark
!!***** STREAM Benchmark license is in MemBenchVector.h file*******!!.
/* Excerpt from STREAM array size notice */
/*
STREAM requires different amounts of memory to run on different
*           systems, depending on both the system cache size(s) and the
*           granularity of the system timer.
*     You should adjust the value of 'STREAM_ARRAY_SIZE' (below)
*           to meet *both* of the following criteria:
*       (a) Each array must be at least 4 times the size of the
*           available cache memory. I don't worry about the difference
*           between 10^6 and 2^20, so in practice the minimum array size
*           is about 3.8 times the cache size.
*           Example 1: One Xeon E3 with 8 MB L3 cache
*               STREAM_ARRAY_SIZE should be >= 4 million, giving
*               an array size of 30.5 MB and a total memory requirement
*               of 91.5 MB.
*           Example 2: Two Xeon E5's with 20 MB L3 cache each (using OpenMP)
*               STREAM_ARRAY_SIZE should be >= 20 million, giving
*               an array size of 153 MB and a total memory requirement
*               of 458 MB.
*       (b) The size should be large enough so that the 'timing calibration'
*           output by the program is at least 20 clock-ticks.
*           Example: most versions of Windows have a 10 millisecond timer
*               granularity.  20 "ticks" at 10 ms/tic is 200 milliseconds.
*               If the chip is capable of 10 GB/s, it moves 2 GB in 200 msec.
*               This means the each array must be at least 1 GB, or 128M elements.
*/

/*********************************************************
              Class implementation
**********************************************************/
mem_bench_impl::AlignMallocBenchmark::AlignMallocBenchmark(_In_ const std::size_t dataSize, _In_ const int nTimes, _In_ const int nThreads) :
m_uidataSize{ dataSize },
m_inTimes{ nTimes },
m_inThreads{ nThreads }
{
	this->m_pdA = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double),32));
	if (this->m_pdA == NULL){
		std::cerr << "Null Pointer: " << std::hex << this->m_pdA << " in Ctor " << std::endl;
		throw std::runtime_error(std::string("Memory allocation failed!! in Ctor: ").append(typeid(AlignMallocBenchmark).name()));
	}
	this->m_pdB = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (this->m_pdB == NULL) {
		std::cerr << "Null Pointer: " << std::hex << this->m_pdB << "in Ctor " << std::endl;
		_mm_free(this->m_pdA);
		throw std::runtime_error(std::string("Memory allocation failed!! in Ctor: ").append(typeid(AlignMallocBenchmark).name()));
	}
	this->m_pdC = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (this->m_pdC == NULL){
		std::cerr << "Null Pointer: " << std::hex << this->m_pdC << "in Ctor: " << std::endl;
		_mm_free(this->m_pdA);
		_mm_free(this->m_pdB);

	}
}

mem_bench_impl::AlignMallocBenchmark::AlignMallocBenchmark(_In_ const AlignMallocBenchmark &rhs) :
m_uidataSize{ rhs.m_uidataSize },
m_inTimes{ rhs.m_inTimes },
m_inThreads{ rhs.m_inThreads }
{
	this->m_pdA = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (this->m_pdA == NULL){
		std::cerr << "Null Pointer: " << std::hex << this->m_pdA << " in Ctor " << std::endl;
		throw std::runtime_error(std::string("Memory allocation failed!! in Ctor: ").append(typeid(AlignMallocBenchmark).name()));
	}
	this->m_pdB = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (this->m_pdB == NULL) {
		std::cerr << "Null Pointer: " << std::hex << this->m_pdB << "in Ctor " << std::endl;
		_mm_free(this->m_pdA);
		throw std::runtime_error(std::string("Memory allocation failed!! in Ctor: ").append(typeid(AlignMallocBenchmark).name()));
	}
	this->m_pdC = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (this->m_pdC == NULL){
		std::cerr << "Null Pointer: " << std::hex << this->m_pdC << "in Ctor: " << std::endl;
		_mm_free(this->m_pdA);
		_mm_free(this->m_pdB);
		throw std::runtime_error(std::string("Memory allocation failed!! in Ctor: ").append(typeid(AlignMallocBenchmark).name()));

	}
	/* Warning !!
	   Assumes that rhs.m_uidataSize == this->m_uidataSize */
	std::copy(&rhs.m_pdA[0], &rhs.m_pdA[0] + rhs.m_uidataSize, &this->m_pdA[0]);
	std::copy(&rhs.m_pdB[0], &rhs.m_pdB[0] + rhs.m_uidataSize, &this->m_pdB[0]);
	std::copy(&rhs.m_pdC[0], &rhs.m_pdC[0] + rhs.m_uidataSize, &this->m_pdC[0]);
}

mem_bench_impl::AlignMallocBenchmark::AlignMallocBenchmark(_In_ AlignMallocBenchmark &&rhs) :
m_uidataSize{ std::move(rhs.m_uidataSize) },
m_inTimes{ std::move(rhs.m_inTimes) },
m_inThreads{ std::move(rhs.m_inThreads) }
{
	this->m_pdA = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (this->m_pdA == NULL){
		std::cerr << "Null Pointer: " << std::hex << this->m_pdA << " in Ctor " << std::endl;
		throw std::runtime_error(std::string("Memory allocation failed!! in Ctor: ").append(typeid(AlignMallocBenchmark).name()));
	}
	this->m_pdB = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (this->m_pdB == NULL) {
		std::cerr << "Null Pointer: " << std::hex << this->m_pdB << "in Ctor " << std::endl;
		_mm_free(this->m_pdA);
		throw std::runtime_error(std::string("Memory allocation failed!! in Ctor: ").append(typeid(AlignMallocBenchmark).name()));
	}
	this->m_pdC = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (this->m_pdC == NULL){
		std::cerr << "Null Pointer: " << std::hex << this->m_pdC << "in Ctor: " << std::endl;
		_mm_free(this->m_pdA);
		_mm_free(this->m_pdB);
		throw std::runtime_error(std::string("Memory allocation failed!! in Ctor: ").append(typeid(AlignMallocBenchmark).name()));

	}
	// Move object state to *this
	// Using std::move which is probably unoptimized for
	// utilization of two load ports.
	std::move(&rhs.m_pdA[0], &rhs.m_pdA[0] + rhs.m_uidataSize, &this->m_pdA[0]);
	std::move(&rhs.m_pdB[0], &rhs.m_pdB[0] + rhs.m_uidataSize, &this->m_pdB[0]);
	std::move(&rhs.m_pdC[0], &rhs.m_pdC[0] + rhs.m_uidataSize, &this->m_pdC[0]);
	// Destroy rhs state.
	_mm_free(rhs.m_pdA); rhs.m_pdA = NULL;
	_mm_free(rhs.m_pdB); rhs.m_pdB = NULL;
	_mm_free(rhs.m_pdC); rhs.m_pdC = NULL;
}

mem_bench_impl::AlignMallocBenchmark::~AlignMallocBenchmark() {

	if (this->m_pdA) _mm_free(this->m_pdA); this->m_pdA = NULL;
	if (this->m_pdB) _mm_free(this->m_pdB); this->m_pdB = NULL;
	if (this->m_pdC) _mm_free(this->m_pdC); this->m_pdC = NULL;
}

/****************************************************************************
      Implementation of AlignMallocBnechmark::run_benchmark()
	                  !!Important!!
	  Calls exit(EXIT_FAILURE) when _mm_malloc fails
*****************************************************************************/
void    mem_bench_impl::AlignMallocBenchmark::run_benchmark() {

	int quantum, BytesPerWord, k;
	double scalar, t, times[4][10];
	double avgtime[4] = { 0.0 };
	const std::string ops_names[4] = { std::string("Copy"), std::string("Scale"),
		std::string("Add"), std::string("Mul-Add") };
	/* Using FLT_MAX because of names conflict when std::numeric_limits<>::max()
	is used.*/
	double mintime[4] = { FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX };
	double maxtime[4] = { 0.0 };
	const double bytes[4] = { 2 * sizeof(double)* static_cast<double>(this->m_uidataSize),
		2 * sizeof(double)* static_cast<double>(this->m_uidataSize),
		3 * sizeof(double)* static_cast<double>(this->m_uidataSize),
		3 * sizeof(double)* static_cast<double>(this->m_uidataSize)
	};
	constexpr int timeSamples{ 20 };
	std::cout << "Memory Benchmark version $Revision: 1.0 $\n" << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	BytesPerWord = sizeof(double);
	std::cout << "This system uses: " << BytesPerWord << " bytes per element." << std::endl;
	std::cout << "In this: " << typeid(AlignMallocBenchmark).name() << " offset is equal: " << 0 << std::endl;
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "operator new allocated: " << " size: " << this->m_uidataSize << " elements with offset: " << 0 << std::endl;
	std::cout << std::fixed << std::setprecision(10) << "Memory per " << typeid(this->m_pdA).name() << ":" << BytesPerWord * (static_cast<double>(this->m_uidataSize) / 1024.0 / 1024.0) << " MiB." << std::endl;
	std::cout << std::fixed << std::setprecision(10) << "Memory per " << typeid(this->m_pdA).name() << ":" << BytesPerWord * (static_cast<double>(this->m_uidataSize) / 1024.0 / 1024.0 / 1024.0) << " GiB." << std::endl;
	std::cout << std::fixed << std::setprecision(10) << "Total memory required: " << (3 * BytesPerWord * (static_cast<double>(this->m_uidataSize) / 1024.0 / 1024.0)) << " MiB." << std::endl;
	std::cout << std::fixed << std::setprecision(10) << "Total memory required: " << (3 * BytesPerWord * (static_cast<double>(this->m_uidataSize) / 1024.0 / 1024.0 / 1024.0)) << " GiB." << std::endl;
	std::cout << "Each computing kernel will be executed: " << this->m_inTimes << " times." << std::endl;
	std::cout << "The *best* time for each kernel (excluding the first iteration) will be used to compute reported BW. \n";
#if (MEM_BENCH_USE_OPENMP) == 0x1
	std::cout << "---------------------------------------------" << std::endl;
	std::cout << "Number of threads requested: " << this->m_inThreads << std::endl;
	k = 0;
#pragma omp parallel
#pragma omp atomic
	k++;
	std::cout << "Number of threads counted: " << k << std::endl;
	/* Intializing Arrays... */
	std::cout << "Initializing  smart pointers to  arrays... ";
	std::size_t dataSize{ this->m_uidataSize };
	const int nThreads{ this->m_inThreads };
	/* Allocate automatic arrays */
	double* a_pdA = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdA == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdA << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::run_benchmark, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		::exit(EXIT_FAILURE);
	}

	double* a_pdB = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdB == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdB << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::run_benchmark, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		_mm_free(a_pdA); a_pdA = NULL;
		::exit(EXIT_FAILURE);
	}

	double* a_pdC = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdC == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdC << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::run_benchmark, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		_mm_free(a_pdA); a_pdA = NULL;
		_mm_free(a_pdB); a_pdB = NULL;
		::exit(EXIT_FAILURE);
	}
	::omp_set_num_threads(nThreads);
#pragma omp parallel for
	for (std::size_t j = 0; j < dataSize; ++j) {
		a_pdA[j] = 1.0;
		a_pdB[j] = 2.0;
		a_pdC[j] = 0.0;
	}
	/* Copy arrays to members */
	std::copy(&a_pdA[0], &a_pdA[0] + dataSize, &this->m_pdA[0]);
	std::copy(&a_pdB[0], &a_pdB[0] + dataSize, &this->m_pdB[0]);
	std::copy(&a_pdC[0], &a_pdC[0] + dataSize, &this->m_pdC[0]);
	std::cout << "Done!!" << std::endl;
	std::cout << "--------------------------------------" << std::endl;
	if ((quantum = mem_bench_impl::MemBenchTimers::check_tick<timeSamples>()) >= 1)
		std::cout << "Your clock granularity appears to be " << quantum << "microsec" << std::endl;
	else {
		std::cout << "Yout clock granularity appears to be less than 1 microsec." << std::endl;
		quantum = 1;
	}
	::omp_set_num_threads(nThreads);
	t = mem_bench_impl::MemBenchTimers::timer_tick();
#pragma omp parallel for
	for (std::size_t j = 0; j < dataSize; ++j)
		a_pdA[j] = 2.0 * a_pdA[j];
	t = 1.0E6 * (mem_bench_impl::MemBenchTimers::timer_tick() - t);
	std::cout << "Each test below will take on the order of: " << static_cast<int>(t) << " microsec." << std::endl;
	std::cout << "Clock ticks: " << static_cast<int>(t / quantum) << std::endl;
	std::cout << "Increase the size of the arrays if this shows that\n";
	std::cout << "you are not getting at least 20 clock ticks per test.\n";

	std::cout << "----------------------------------------------------" << std::endl;
	std::cout << "WARNING -- The above is only a rough guideline.\n";
	std::cout << "For best results, please be sure you know the\n";
	std::cout << "precision of your system timer.\n";
	std::cout << "----------------------------------------------------" << std::endl;
	std::cout << "WARNING -- The above is only a rough guideline.\n";
	std::cout << "For best results, please be sure you know the\n";
	std::cout << "precision of your system timer.\n";

	std::cout << "----------------------------------------------------" << std::endl;

	/*	--- MAIN LOOP --- repeat test cases NTIMES times --- */
	scalar = 3.0;
	for (int k{ 0 }; k != this->m_inTimes; ++k) {

		/* Calling copy_operation NTIMES */
		times[0][k] = mem_bench_impl::MemBenchTimers::timer_tick();
		copy_operation();
		times[0][k] = mem_bench_impl::MemBenchTimers::timer_tick() - times[0][k];

		/* Calling scaling operation NTIMES */
		times[1][k] = mem_bench_impl::MemBenchTimers::timer_tick();
		scaling_operation(scalar);
		times[1][k] = mem_bench_impl::MemBenchTimers::timer_tick() - times[1][k];

		/* Calling addition operations NTIMES */
		times[2][k] = mem_bench_impl::MemBenchTimers::timer_tick();
		addition_operation();
		times[2][k] = mem_bench_impl::MemBenchTimers::timer_tick() - times[2][k];

		/* Calling mul_add_operations NTIMES */
		times[3][k] = mem_bench_impl::MemBenchTimers::timer_tick();
		mul_add_operation(scalar);
		times[3][k] = mem_bench_impl::MemBenchTimers::timer_tick() - times[3][k];
	}
	/* Summary */
	for (int k = 1; k != this->m_inTimes; ++k){
		for (int j = 0; j != 4; ++j) {
			avgtime[j] = avgtime[j] + times[j][k];
			mintime[j] = min(mintime[j], times[j][k]);
			maxtime[j] = max(maxtime[j], times[j][k]);
		}
	}

	std::cout << "  Operation        Best Rate MB/s         Avg time        Min time          Max time" << std::endl;
	for (int j{ 0 }; j < 4; ++j) {
		avgtime[j] = avgtime[j] / static_cast<double>(this->m_inTimes - 1);

		std::cout << std::setw(2) << ops_names[j].c_str() << std::setw(25) << 1.0E-6*bytes[j] / mintime[j] << std::setw(19) << avgtime[j] << std::setw(21) << mintime[j] << std::setw(24) << maxtime[j] << std::endl;

	}
	std::cout << "--------------------------------------------" << std::endl;

	/* Check results */
	check_bench_results();
	/* Deallocating automatic arrays */
	_mm_free(a_pdA);
	_mm_free(a_pdB);
	_mm_free(a_pdC);
#else
                             /* Using single thread only */
    std::cout << "Using single Thread -- Main Thread of ID:" << ::GetCurrentThreadId() << std::endl;
   for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j) {
	   this->m_pdA[j] = 1.0;
	   this->m_pdB[j] = 2.0;
	   this->m_pdC[j] = 0.0;
   }

std::cout << "--------------------------------------------" << std::endl;

   if ((quantum = mem_bench_impl::MemBenchTimers::check_tick<timeSamples>()) >= 1)
        std::cout << "Your clock granularity/precision appears to be: " << quantum << std::endl;
   else {
	   std::cout << "Your clock granularity appears to be less than one microsec" << std::endl;
	   quantum = 1;
    }

   t = mem_bench_impl::MemBenchTimers::timer_tick();
   for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j)
        this->m_pdA[j] = 2.0 * this->m_pdA[j];
   t = 1.0E6 * (mem_bench_impl::MemBenchTimers::timer_tick() - t);

std::cout << "Each test below will take on the order of: " << static_cast<int>(t) << " microsec." << std::endl;
std::cout << "Clock ticks: " << static_cast<int>(t / quantum) << std::endl;
std::cout << "Increase the size of the arrays if this shows that\n";
std::cout << "you are not getting at least 20 clock ticks per test.\n";

std::cout << "--------------------------------------------" << std::endl;
std::cout << "WARNING -- The above is only a rough guideline.\n";
std::cout << "For best results, please be sure you know the\n";
std::cout << "precision of your system timer.\n";
std::cout << "--------------------------------------------" << std::endl;

/*	--- MAIN LOOP --- repeat test cases NTIMES times --- */
      scalar = 3.0;
for (int k{ 0 }; k != this->m_inTimes; ++k) {

	//std::cout << " copy_operation NTIMES " << std::endl;
	times[0][k] = mem_bench_impl::MemBenchTimers::timer_tick();
	copy_operation();
	times[0][k] = mem_bench_impl::MemBenchTimers::timer_tick() - times[0][k];

	/* Calling scaling operation NTIMES */
	times[1][k] = mem_bench_impl::MemBenchTimers::timer_tick();
	scaling_operation(scalar);
	times[1][k] = mem_bench_impl::MemBenchTimers::timer_tick() - times[1][k];

	/* Calling addition operations NTIMES */
	times[2][k] = mem_bench_impl::MemBenchTimers::timer_tick();
	addition_operation();
	times[2][k] = mem_bench_impl::MemBenchTimers::timer_tick() - times[2][k];

	/* Calling mul_add_operations NTIMES */
	times[3][k] = mem_bench_impl::MemBenchTimers::timer_tick();
	mul_add_operation(scalar);
	times[3][k] = mem_bench_impl::MemBenchTimers::timer_tick() - times[3][k];
}

for (int k = 1; k != this->m_inTimes; ++k){
	for (int j = 0; j != 4; ++j) {
		avgtime[j] = avgtime[j] + times[j][k];
		mintime[j] = min(mintime[j], times[j][k]);
		maxtime[j] = max(maxtime[j], times[j][k]);
	}
}

std::cout << "  Operation        Best Rate MB/s         Avg time        Min time          Max time" << std::endl;
for (int j{ 0 }; j < 4; ++j) {
	avgtime[j] = avgtime[j] / static_cast<double>(this->m_inTimes - 1);

	std::cout << std::setw(2) << ops_names[j].c_str() << std::setw(25) << 1.0E-6*bytes[j] / mintime[j] << std::setw(19) << avgtime[j] << std::setw(21) << mintime[j] << std::setw(24) << maxtime[j] << std::endl;

}
std::cout << "--------------------------------------------" << std::endl;

/* Check results */
check_bench_results();

#endif

}


/*****************************************************************************
       Implementation of AlignMallocBenchmar::check_bench_results()
******************************************************************************/
void  mem_bench_impl::AlignMallocBenchmark::check_bench_results() {

	/* Reproducing initialization values */
	double aj{ 1.0 }, bj{ 2.0 }, cj{ 0.0 }, scalar{ 3.0 };
	double aSumErr{ 0.0 }, bSumErr{ 0.0 }, cSumErr{ 0.0 };
	double aAvgErr, bAvgErr, cAvgErr;
	int err{ 0 }, ierr;
	constexpr double EPS{ 1.0E-13 };
	aj = 2.0E0 * aj;
#ifndef MEM_BENCH_VERBOSE
	std::cout << "Entering: UniquePtrBenchmark::check_bench_results() " << std::endl;
#endif
	for (int k{ 0 }; k != this->m_inTimes; ++k){

		cj = aj;
		bj = scalar * cj;
		cj = aj + bj;
		aj = bj + scalar*cj;
	}

	/* accumulate deltas between observed and expected results */
	for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j) {

		aSumErr += std::abs(this->m_pdA[j] - aj);
		bSumErr += std::abs(this->m_pdB[j] - bj);
		cSumErr += std::abs(this->m_pdC[j] - cj);
	}
	aAvgErr = aSumErr / static_cast<double>(this->m_uidataSize);
	bAvgErr = bSumErr / static_cast<double>(this->m_uidataSize);
	cAvgErr = cSumErr / static_cast<double>(this->m_uidataSize);
	/* Epsilon is set to 1.0E-13 */
	if (std::abs(aAvgErr / aj) > EPS){
		err++;
		std::cout << "Failed Validation on malloc'ed m_pdA, AvgRelAbsErr > EPS: " << std::setprecision(13) << std::fixed << EPS << std::endl;
		std::cout << "Expected value: " << std::fixed << std::showpoint << std::setprecision(13) <<
			aj << std::endl;
		std::cout << "AvgAbsErr: " << std::fixed << std::showpoint << std::setprecision(13) << aAvgErr << std::endl;
		std::cout << "AvgRelAbsErr: " << std::fixed << std::showpoint << std::setprecision(13) << std::abs(aAvgErr / aj) << std::endl;
		ierr = 0;
		for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j) {
			if (std::abs(this->m_pdA[j] / aj - 1.0) > EPS) {
				ierr++;
#ifdef MEM_BENCH_VERBOSE
				if (ierr < 10){
					std::cout << std::fixed << std::showpoint << std::setprecision(13) << "std::unique_ptr m_spA: index " << j << " expected: " << aj << " observed: " << this->m_pdA[j] <<
						"rel error: " << std::abs((aj - this->m_pdA[j]) / aAvgErr) << std::endl;
				}
#endif
			}
		}
		std::cout << "For malloc'ed m_pdA: " << ierr << " errors were found." << std::endl;
	}

	if (std::abs(bAvgErr / bj) > EPS) {
		err++;
		std::cout << "Failed Validation on malloc'ed m_pdB, AvgRelAbsErr > EPS: " << std::setprecision(13) << std::fixed << EPS << std::endl;
		std::cout << "Expected value: " << std::fixed << std::showpoint << std::setprecision(13) <<
			bj << std::endl;
		std::cout << "AvgAbsErr: " << std::fixed << std::showpoint << std::setprecision(13) << bAvgErr << std::endl;
		std::cout << "AvgRelAbsErr: " << std::fixed << std::showpoint << std::setprecision(13) << std::abs(bAvgErr / bj) << std::endl;
		ierr = 0;
		for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j) {
			if (std::abs(this->m_pdB[j] / bj - 1.0) > EPS) {
				ierr++;
#ifdef MEM_BENCH_VERBOSE
				if (ierr < 10) {
					std::cout << std::fixed << std::showpoint << std::setprecision(13) << "std::unique_ptr m_spB: index " << j << " expected: " << bj << "observed: " << this->m_pdB[j] <<
						"rel error: " << std::abs((bj - this->m_pdB[j]) / bAvgErr) << std::endl;
				}
#endif
			}
		}
		std::cout << "For malloc'ed m_pdB: " << ierr << " errors were found." << std::endl;
	}

	if (std::abs(cAvgErr / cj) > EPS) {
		err++;
		std::cout << "Failed Validation on malloc'ed m_pdC, AvgRelAbsErr > EPS: " << std::setprecision(13) << std::fixed << EPS << std::endl;
		std::cout << "Expected value: " << std::fixed << std::showpoint << std::setprecision(13) <<
			cj << std::endl;
		std::cout << "AvgAbsErr: " << std::fixed << std::showpoint << std::setprecision(13) << cAvgErr << std::endl;
		std::cout << "AvgRelAbsErr: " << std::fixed << std::showpoint << std::setprecision(13) << std::abs(cAvgErr / cj) << std::endl;
		ierr = 0;
		for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j) {
			if (std::abs(this->m_pdC[j] / cj - 1.0) > EPS) {
				ierr++;
#ifdef MEM_BENCH_VERBOSE
				if (ierr < 10) {
					std::cout << std::fixed << std::showpoint << std::setprecision(13) << " std::unique_ptr m_aC: index " << j << " expected: " << cj << "observed: " << this->m_pdC[j] <<
						"rel error: " << std::abs((cj - this->m_pdC[j]) / cAvgErr) << std::endl;
				}
#endif
			}
		}
		std::cout << "For malloc'ed m_pdC: " << ierr << " errors were found." << std::endl;

	}
	if (err == 0) {
		std::cout << std::fixed << std::showpoint << std::setprecision(13) << "Solution Validates: average error less than: " << EPS << " on all 3 malloc'ed arrays. " << std::endl;
	}
#ifndef MEM_BENCH_DEBUG_VERBOSE
	std::cout << "Exiting: AlignMallocBenchmark::check_bench_results() " << std::endl;
#endif
}

/******************************************************************************
          Implementation of AlignMallocBenchmark::copy_operation()
		                   !!Important!!
		  Calls exit(EXIT_FAILURE) when _mm_malloc fails
*******************************************************************************/
void     mem_bench_impl::AlignMallocBenchmark::copy_operation() {
	/*Uses OpenMP this will imply creation of automatic
	variables and additional copy operations*/
#if (MEM_BENCH_USE_OPENMP) == 0x1
	/* Warning possible overhead of 3 vector copy-ctors
	Accurate time of these copy operations shall
	be measured and subtracted from the results*/
	/* Handle memory allocation errors */
	double*  a_pdA = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdA == NULL){
		std::cerr << "Null Pointer: " << std::hex << a_pdA << " at loc #: " << __LINE__ << std::endl;
		std::cerr << " in AlignMallocBenchmark::copy_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		::exit(EXIT_FAILURE);
	}

	double* a_pdC = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdC == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdC << " at loc #: " << __LINE__ << std::endl;
		std::cerr << " in AlignMallocBenchmark::copy_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		_mm_free(a_pdA); a_pdA = NULL;

		::exit(EXIT_FAILURE);
	}

	
	/* copy object state to automatic arrays. */
	/* Using unoptimized scalar version , moving one double/cycle/port with vmovsd */
	std::copy(&this->m_pdA[0], &this->m_pdA[0] + this->m_uidataSize, &a_pdA[0]);
	std::copy(&this->m_pdC[0], &this->m_pdC[0] + this->m_uidataSize, &a_pdC[0]);
	const std::size_t dataSize{ this->m_uidataSize };
	const int nThreads{ this->m_inThreads };
	::omp_set_num_threads(nThreads);
#pragma omp parallel for
	for (std::size_t j = 0; j < dataSize; ++j)
		a_pdC[j] = a_pdA[j];
	
	std::copy(&a_pdC[0], &a_pdC[0] + dataSize, &this->m_pdC[0]);
	// Deallocate automatic arrays!!
	_mm_free(a_pdA); 
	_mm_free(a_pdC); 
#else
	/* Uses single threaded version. Run this first
	take timing and compare it with multi-threaded
	version*/
	for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j)
		this->m_pdC[j] = this->m_pdA[j];
#endif
#ifndef MEM_BENCH_DEBUG_VERBOSE
	std::cout << "Executing in: UniquePtrBenchmark::copy_operation()" << std::endl;
#endif
}

/******************************************************************************
        Implementation of AlignMallocBenchmark::scaling_operation()
		                    !!Important!!
		Calls exit(EXIT_FAILURE) when _mm_malloc fails
*******************************************************************************/
void     mem_bench_impl::AlignMallocBenchmark::scaling_operation(_In_ const double scalar) {
	/*Uses OpenMP this will imply creation of automatic
	variables and additional copy operations*/
#if (MEM_BENCH_USE_OPENMP) == 0x1
	/* Warning possible overhead of 3 vector copy-ctors
	Accurate time of these copy operations shall
	be measured and subtracted from the results*/
	/* Handle memory allocation errors */
	double* a_pdB = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdB == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdB << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBenchmark::scaling_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		::exit(EXIT_FAILURE);
	}

	double* a_pdC = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdC == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdC << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::scaling_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		_mm_free(a_pdB); a_pdB = NULL;
		::exit(EXIT_FAILURE);
	}
	/* copy object state to automatic arrays. */
	/* Using unoptimized scalar version , moving one double/cycle/port with vmovsd */
	std::copy(&this->m_pdB[0], &this->m_pdB[0] + this->m_uidataSize, &a_pdB[0]);
	std::copy(&this->m_pdC[0], &this->m_pdC[0] + this->m_uidataSize, &a_pdC[0]);
	const std::size_t dataSize{ this->m_uidataSize };
	const int nThreads{ this->m_inThreads };
	::omp_set_num_threads(nThreads);
#pragma omp parallel for
	for (std::size_t j = 0; j < dataSize; ++j)
		a_pdB[j] = scalar * a_pdC[j];

	std::copy(&a_pdB[0], &a_pdB[0] + dataSize, &this->m_pdB[0]);
	// Deallocate automatic arrays!!
	_mm_free(a_pdB); 
	_mm_free(a_pdC);
#else
	/* Uses single threaded version. Run this first
	take timing and compare it with multi-threaded
	version*/
	for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j)
		this->m_pdB[j] = scalar * this->m_pdC[j];
#endif
#ifndef MEM_BENCH_DEBUG_VERBOSE
	std::cout << "Executing in: AlignMallocBenchmark::scaling_operation()" << std::endl;
#endif
}

/*******************************************************************************
      Implementation of AlignMalocBenchmark::addition_operation()
	                      !!Important!!
		  Calls exit(EXIT_FAILURE) when _mm_malloc fails
********************************************************************************/
void   mem_bench_impl::AlignMallocBenchmark::addition_operation() {
	/*Uses OpenMP this will imply creation of automatic
	variables and additional copy operations*/
#if (MEM_BENCH_USE_OPENMP) == 0x1
	/* Warning possible overhead of 3 vector copy-ctors
	Accurate time of these copy operations shall
	be measured and subtracted from the results*/
	/* Handle memory allocation errors */
	double* a_pdA = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdA == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdA << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::addition_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		::exit(EXIT_FAILURE);
	}

	double* a_pdB = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdB == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdB << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::addition_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		_mm_free(a_pdA); a_pdA = NULL;
		::exit(EXIT_FAILURE);
	}

	double* a_pdC = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdC == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdC << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::addition_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		_mm_free(a_pdA); a_pdA = NULL;
		_mm_free(a_pdB); a_pdB = NULL;
		::exit(EXIT_FAILURE);
	}
	/* copy object state to automatic arrays. */
	/* Using unoptimized scalar version , moving one double/cycle/port with vmovsd */
	std::copy(&this->m_pdA[0], &this->m_pdA[0] + this->m_uidataSize, &a_pdA[0]);
	std::copy(&this->m_pdB[0], &this->m_pdB[0] + this->m_uidataSize, &a_pdB[0]);
	std::copy(&this->m_pdC[0], &this->m_pdC[0] + this->m_uidataSize, &a_pdC[0]);
	const std::size_t dataSize{ this->m_uidataSize };
	const int nThreads{ this->m_inThreads };
	::omp_set_num_threads(nThreads);
#pragma omp parallel for
	for (std::size_t j = 0; j < dataSize; ++j)
		a_pdC[j] = a_pdA[j] + a_pdB[j];

	std::copy(&a_pdC[0], &a_pdC[0] + dataSize, &this->m_pdC[0]);
	// Deallocate automatic arrays !!
	_mm_free(a_pdA);
	_mm_free(a_pdB);
	_mm_free(a_pdC);
#else
	/* Uses single threaded version. Run this first
	take timing and compare it with multi-threaded
	version*/
	for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j)
		this->m_pdC[j] = this->m_pdA[j] + this->m_pdB[j];
#endif
#ifndef MEM_BENCH_DEBUG_VERBOSE
	std::cout << "Executing in: AlignMallocBenchmark::addition_operation()" << std::endl;
#endif
}

/***************************************************************************
      Implementation of AlignMallocBnechmark::mul_add_operation()
	                    !!Important!!
		  Calls exit(EXIT_FAILURE) when _mm_malloc fails
****************************************************************************/
void    mem_bench_impl::AlignMallocBenchmark::mul_add_operation(_In_ const double scalar) {
	/*Uses OpenMP this will imply creation of automatic
	variables and additional copy operations*/
#if (MEM_BENCH_USE_OPENMP) == 0x1
	/* Warning possible overhead of 3 vector copy-ctors
	Accurate time of these copy operations shall
	be measured and subtracted from the results*/
	/* Handle memory allocation errors */
	double* a_pdA = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdA == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdA << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::mul_add_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		::exit(EXIT_FAILURE);
	}

	double* a_pdB = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdB == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdB << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::mul_add_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		_mm_free(a_pdA); a_pdA = NULL;
		::exit(EXIT_FAILURE);
	}

	double* a_pdC = reinterpret_cast<double*>(_mm_malloc(this->m_uidataSize * sizeof(double), 32));
	if (a_pdC == NULL) {
		std::cerr << "Null Pointer: " << std::hex << a_pdC << " at loc #: " << __LINE__;
		std::cerr << " in AlignMallocBnechmark::mul_add_operation, cleaning state and calling exit(EXIT_FAILURE)" << std::endl;
		_mm_free(a_pdA); a_pdA = NULL;
		_mm_free(a_pdB); a_pdB = NULL;
		::exit(EXIT_FAILURE);
	}
	/* copy object state to automatic arrays. */
	/* Using unoptimized scalar version , moving one double/cycle/port with vmovsd */
	std::copy(&this->m_pdA[0], &this->m_pdA[0] + this->m_uidataSize, &a_pdA[0]);
	std::copy(&this->m_pdB[0], &this->m_pdB[0] + this->m_uidataSize, &a_pdB[0]);
	std::copy(&this->m_pdC[0], &this->m_pdC[0] + this->m_uidataSize, &a_pdC[0]);
	const std::size_t dataSize{ this->m_uidataSize };
	const int nThreads{ this->m_inThreads };
	::omp_set_num_threads(nThreads);
#pragma omp parallel for
	for (std::size_t j = 0; j < dataSize; ++j)
		a_pdA[j] = a_pdB[j] + scalar * a_pdC[j];

	std::copy(&a_pdA[0], a_pdA[0] + dataSize, &this->m_pdA[0]);
	// Deallocate automatic arrays !!
	_mm_free(a_pdA);
	_mm_free(a_pdB);
	_mm_free(a_pdC);
#else
	/* Uses single threaded version. Run this first
	take timing and compare it with multi-threaded
	version*/
	for (std::size_t j{ 0 }; j != this->m_uidataSize; ++j)
		this->m_pdA[j] = this->m_pdB[j] + scalar * this->m_pdC[j];
#endif
#ifndef MEM_BENCH_DEBUG_VERBOSE
	std::cout << "Executing in: AlignMallocBenchmark::mul_add_operation()" << std::endl;
#endif
}