
#ifndef __CUDA_TIME_WAIT_H__
#define __CUDA_TIME_WAIT_H__

#include <time.h>

void wait (int milliseconds) {
  clock_t endwait;
  const int milli = 1000;
  endwait = clock() + milliseconds * CLOCKS_PER_SEC / milli;
  while (clock() < endwait) {}
}

#endif // __CUDA_TIME_WAIT_H__
