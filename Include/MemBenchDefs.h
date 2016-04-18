#ifndef _MEM_BENCH_DEFS_H_04_16_16
#define _MEM_BENCH_DEFS_H_04_16_16 0x1

/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Memory Benchmarking tool adapted from John D. McCalpin "STREAM" benchmark
@file MemBenchDefs.h
@author: Bernard Gingold
@version:  1.0.0  15/04/2016 19:59
@description: MemBenchDefs.h
@reference: none

  *****************NOTICE!!***********************
  STREAM License is provided by the implementation
                  header files.
  *************************************************
*/

/************************************************
   Frequently used header files.
*************************************************/

#include <vector>
#include <valarray>
#include <array>
#include <memory>
#include <iostream>
#include <iomanip>
//#include <processthreadsapi.h>
#include <Windows.h>

#if defined __INTEL_COMPILER
#include <omp.h>
#else
#include <omp.h>
#endif

#ifndef MEM_BENCH_USE_OPENMP
#define MEM_BENCH_USE_OPENMP 0x0
#endif
#ifndef MEM_BENCH_VERBOSE
#define MEM_BENCH_VERBOSE 0x1
#endif
#ifndef MEM_BENCH_DEBUG_VERBOSE
#define MEM_BENCH_DEBUG_VERBOSE 0x1
#endif


#endif /*_MEM_BENCH_DEFS_H_04_16_16*/