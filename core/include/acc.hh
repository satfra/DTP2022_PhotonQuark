#pragma once

// OpenACC decoration helper.
//
// ACC_ROUTINE_SEQ expands to `#pragma acc routine seq` only when compiling with
// an OpenACC-enabled compiler (nvc++ -acc defines the standard _OPENACC macro).
// Under g++/clang (CPU/OpenMP build) it expands to nothing, so -Wall's
// -Wunknown-pragmas stays quiet and the header is a no-op. Place it immediately
// before a function (or member function) definition to make that function
// callable from within an OpenACC compute region.
#ifdef _OPENACC
  #define ACC_ROUTINE_SEQ _Pragma("acc routine seq")
#else
  #define ACC_ROUTINE_SEQ
#endif
