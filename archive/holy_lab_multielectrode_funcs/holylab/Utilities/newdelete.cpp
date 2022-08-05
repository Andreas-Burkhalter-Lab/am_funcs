//
// new and delete C++ operators for MATLAB.
//
// Compile this file with your C++ MEX file
//
// Petter Strandmark 2009
//

#include "newdelete.h"

#include "mex.h"


	
void *operator new(size_t size)
{
	void* ptr = mxMalloc(size);
	if (!ptr) {
		throw std::bad_alloc();
	}
}
void *operator new(size_t size, const std::nothrow_t&) throw()
{
	return mxMalloc(size);
}
void *operator new[](size_t size)
{
	void* ptr = mxMalloc(size);
	if (!ptr) {
		throw std::bad_alloc();
	}
}
void *operator new[](size_t size, const std::nothrow_t&) throw()
{
	return mxMalloc(size);
}

void operator delete(void *p)
{
	if (p) {
		return mxFree(p);
	}
}
void operator delete[](void *p)
{
	if (p) {
		return mxFree(p);
	}
}
