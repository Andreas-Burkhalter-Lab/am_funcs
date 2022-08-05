//
// new and delete C++ operators for MATLAB.
//
// Include this file in your C++ MEX file
//
// Petter Strandmark 2009
//

#ifndef NEWDELETE_HEADER
#define NEWDELETE_HEADER

#include <new>



void* operator new(size_t size);
void* operator new(size_t size, const std::nothrow_t&) throw();
void* operator new[](size_t size);
void* operator new[](size_t size, const std::nothrow_t&) throw();

void operator delete(void *p);
void operator delete[](void *p);

#endif
