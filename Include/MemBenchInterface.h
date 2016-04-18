#ifndef _MEM_BENCH_INTERFACE_H_04_15_16
#define _MEM_BENCH_INTERFACE_H_04_15_16

/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Memory Benchmarking tool adapted from John D. McCalpin "STREAM" benchmark
@file MemBenchInterface.h.h
@author: Bernard Gingold
@version:  1.0.0  15/04/2016 19:59
@description: MemBenchInterface.h
@reference: STREAM benchmark.

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

#include "MemBenchDefs.h"


namespace mem_bench_iface {

	/********************************************
	   Abstarct class which acts as an interface 
	   to classes which implements pure abstract
	   functions.
	   Classes are based on STL data containers,
	   naked new allocation, smart pointer wrapped
	   dynamic allocation and C-like allocation
	   based on calls to _mm_malloc
	*********************************************/
	class MemBenchIface {

	public:

		/*******************************
		  Pure abstract functions 
		  overriden by the subclasses
		********************************/

		/* Pure abstract functions. Runs memory benchmark
		   overriden by specific subclass .*/

		virtual void   run_benchmark() = 0;

		/* Pure abstract function. Checks benchmark results.
		   Overriden by the subclass implementation.*/

		virtual void   check_bench_results() = 0;

		/* Pure abstract function. Runs simple computing
		   kernel operation: copying.
		   Overriden by the subclass implementation.*/
		virtual void    copy_operation() = 0;

		/* Pure abstract function. Runs simple computing
		   kernel operation: vector scaling by scalar.
		   Overriden by the subclass implementation.*/
		virtual void    scaling_operation(_In_ const double) = 0;

		/* Pure abstract function. Runs simple computing
		   kernel operation: vector addition.
		   Overriden by the subclass implementation.*/
		virtual void    addition_operation() = 0;

		/* Pure abstract function. Runs simple computing
		   kernel operation: vector-scalar mul-add.
		   Overriden by the subclass implementation. */
		virtual void    mul_add_operation(_In_ const double) = 0;

		/* Vurtual Destructor*/
		virtual  ~MemBenchIface() {}
			
	};


	
}
#endif /*_MEM_BENCH_INTERFACE_H_04_15_16*/