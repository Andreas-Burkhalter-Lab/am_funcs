#include <time.h>

#include "timer_g.h"

#ifdef _UNIX_G
#include <sys/times.h>
#include <sys/param.h>

struct tms startpoint, readpoint;
void start_timer()
{
   times(&startpoint);
}

/* /------------------------------------ */
float read_timer()
{
   clock_t elapsed;

   times(&readpoint);
   elapsed=readpoint.tms_utime + readpoint.tms_cutime
          -startpoint.tms_utime- startpoint.tms_cutime;

   return ((float)elapsed)/HZ;//CLOCKS_PER_SEC;

//   return ((float)elapsed)/CLOCKS_PER_SEC;

}

#elif defined(_MULTITHREADED)
#include <sys/time.h>

struct timeval startpoint, readpoint;
struct timezone ttz;
void start_timer()
{
   gettimeofday(&startpoint, &ttz);
}

/* /------------------------------------ */
float read_timer()
{
   gettimeofday(&readpoint, &ttz);

   return (readpoint.tv_sec - startpoint.tv_sec)
          +(readpoint.tv_usec - startpoint.tv_usec)/1.0e6;

}

#else

clock_t startpoint, readpoint;
void start_timer()
{
   startpoint=clock();
}

float read_timer()
{
   clock_t elapsed;

   readpoint=clock();
   elapsed=readpoint-startpoint;

   return ((float)elapsed)/CLOCKS_PER_SEC;

}

#endif


