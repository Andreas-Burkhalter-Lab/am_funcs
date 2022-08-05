#ifndef timer_g_h
#define timer_g_h

#ifdef __cplusplus
   extern "C" {
#endif


/* todo: temp use multithread */
#define _MULTITHREADED

void start_timer();
float read_timer();/*return sec*/


#ifdef __cplusplus
   }
#endif

#endif
