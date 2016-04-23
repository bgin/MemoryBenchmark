#ifndef _MEM_BENCH_ALIGN_MALLOC_H_04_22_16
#define _MEM_BENCH_ALIGN_MALLOC_H_04_22_16

/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Memory Benchmarking tool adapted from John D. McCalpin "STREAM" benchmark.
@file MemBenchAlignMalloc.h
@author: Bernard Gingold
@version:  1.0 15/04/2016 19:59
@description: MemBenchAlignMalloc.h
@reference: see below

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

#include "MemBenchInterface.h"
#include "MemBenchHelpers.h"

namespace mem_bench_impl {

	/***************************************************
	   Memory Benchmark implementation based on calls
	   to _mm_malloc in order to allocate the memory
	   aligned on 32-byte boundary.
	****************************************************/

	/********************************************
	            Warning!!
	      Derived class is final
	********************************************/
	class AlignMallocBenchmark final : mem_bench_iface::MemBenchIface {


		/******************************************
		       Constructors and Destructors
		******************************************/

	public:

		/* Explicitly defult Ctor - do not use it
		   object in non/null initialized state*/
		AlignMallocBenchmark() = default;

		/* "Main" class Ctor */
		AlignMallocBenchmark(_In_ const std::size_t, _In_ const int, _In_ const int) noexcept(false);

		/* Copy -  Ctor*/
		AlignMallocBenchmark(_In_ const AlignMallocBenchmark &)noexcept(false);

		/* Move - Ctor */
		AlignMallocBenchmark(_In_ AlignMallocBenchmark &&rhs)noexcept(false);

		/* Class Dtor */
		~AlignMallocBenchmark();
		
		/* Pure abstract functions. Runs memory benchmark
		overriden here .
		*/
		virtual    void   run_benchmark() override;

		/* Pure abstract function. Checks benchmark results.
		Overriden here.*/
		virtual    void   check_bench_results() override;

		/* Pure abstract function. Runs simple computing
		kernel operation: copying.
		Overriden here.*/
		virtual    void   copy_operation() override;

		/* Pure abstract function. Runs simple computing
		kernel operation: vector scaling by scalar.
		Overriden here.*/
		virtual    void   scaling_operation(_In_ const double) override;

		/* Pure abstract function. Runs simple computing
		kernel operation: vector addition.
		Overriden here.*/
		virtual    void    addition_operation() override;

		/* Pure abstract function. Runs simple computing
		kernel operation: vector-scalar mul-add.
		Overriden here. */
		virtual    void   mul_add_operation(_In_ const double) override;



	private:

		/* class member m_uidataSize denotes amount of reserved memory/elements*/
		std::size_t  m_uidataSize;

		/* class member m_inTimes denotes number of benchmark runs.*/
		int  m_inTimes;

		/* class member m_inTimes denotes number of OpenMP threads scheduled
		to run computation kernels.*/
		int  m_inThreads;

		/* class member m_pdX denotes dynamically allocated array
		   by the call to _mm_malloc */
		double*  m_pdA;

		/*
		 class member m_pdB denotes dynamically allocated array
		   by the call to _mm_malloc. */
		double*  m_pdB;

		/*
		class member m_pdC denotes dynamically allocated array
		by the call to _mm_malloc. */
		double*  m_pdC;

	};
}
#endif /*_MEM_BENCH_ALIGN_MALLOC_H_04_22_16*/