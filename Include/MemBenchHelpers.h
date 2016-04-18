#ifndef _MEM_BENCH_HELPERS_H_04_15_16
#define _MEM_BENCH_HELPERS_H_04_15_16 0x1

/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Memory Benchmarking tool adapted from John D. McCalpin "STREAM" benchmark
@file MemBenchHelpers.h.h
@author: Bernard Gingold
@version:  1.0.0  15/04/2016 19:59
@description: MemBenchHelpers.h
@reference: none

STREAM License is provided here:
License:                                                              */
/*  1. You are free to use this program and/or to redistribute           */
/*     this program.                                                     */
/*  2. You are free to modify this program for your own use,             */
/*     including commercial use, subject to the publication              */
/*     restrictions in item 3.                                           */
/*  3. You are free to publish results obtained from running this        */
/*     program, or from works that you derive from this program,         */
/*     with the following limitations:                                   */
/*     3a. In order to be referred to as "STREAM benchmark results",     */
/*         published results must be in conformance to the STREAM        */
/*         Run Rules, (briefly reviewed below) published at              */
/*         http://www.cs.virginia.edu/stream/ref.html                    */
/*         and incorporated herein by reference.                         */
/*         As the copyright holder, John McCalpin retains the            */
/*         right to determine conformity with the Run Rules.             */
/*     3b. Results based on modified source code or on runs not in       */
/*         accordance with the STREAM Run Rules must be clearly          */
/*         labelled whenever they are published.  Examples of            */
/*         proper labelling include:                                     */
/*           "tuned STREAM benchmark results"                            */
/*           "based on a variant of the STREAM benchmark code"           */
/*         Other comparable, clear, and reasonable labelling is          */
/*         acceptable.                                                   */
/*     3c. Submission of results to the STREAM benchmark web site        */
/*         is encouraged, but not required.                              */
/*  4. Use of this program or creation of derived works based on this    */
/*     program constitutes acceptance of these licensing restrictions.   */
/*  5. Absolutely no warranty is expressed or implied.


*/

#include <algorithm>
#include <math.h>
#include <Windows.h>/* For QueryPerformanceCounter/Frequency */

namespace mem_bench_impl {



	/********************************************
	    Helper class static functions for 
		time measurement and time measurement
		granularity.
	*********************************************/
	class MemBenchTimers {


	public:

		/*  Collect a sequence of M  unique time values 
		    from the system. */
		template<const int M> static int check_tick();

		/* Wall clock timer value for Windows systems
		 only. Removed include of <sys/time.h> 
		 Originally this function was named "mysecond()"*/
		static  auto   timer_tick()->double;
		
		
	};

#include "MemBenchHelpers.inl"

}
#endif /*_MEM_BENCH_HELPERS_H_04_15_16*/