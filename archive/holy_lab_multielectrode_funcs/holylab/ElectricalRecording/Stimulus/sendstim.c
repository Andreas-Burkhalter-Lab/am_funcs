/* 
 * sendstim: send a sequence of commands to the valvecontroller
 *
 * This reads a valve sequence from its standard input and then
 * schedules a series of valve switches.  The input has the following form:
 *   n
 *   valve1 time1
 *   valve2 time2
 *   ...
 *   valven timen
 * where "n" on the first line is the number of lines that follow.
 *
 * Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>
 */

#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include <comedilib.h>
#include <stdlib.h>
#include <sched.h>
#include <sys/mman.h>

#define MAXLINE 4096

/* define COMEDI to 1 if you want to generate the DIO output,
   and to 0 if you just want to test timing */
#define COMEDI 1
/* define REALTIME to 1 if you want to try to run this in realtime.
   You'll need root privileges.  In my experience, this has made no 
   difference.*/
#define REALTIME 0

void getcurrenttime(struct timeval *tvP,char *errmsg);
int set_realtime_priority(void);

int main(int argc, char *argv[])
{
  char     line[MAXLINE];
  int      *valveseq,ntransitions;
  float    *timeseq;
  int      valvenum,savevalvenum;
  float    timeswitch;
  double   ctimeoffset,sleeptime;
  struct timeval starttime,curtime;
  /* These variables are needed to control the DIO lines on the
     National Instruments PCI-6071E. */
  char     *comdevicename = "/dev/comedi0";
  unsigned int subdevice = 2;
  int ndiochannels = 4;
  unsigned int writemask;
  unsigned int bits;
  comedi_t *dev;
  int      stype;
  int      scanret;
  int      i;

  /* Since all times are referenced to 0 = starttime, figure this
     out right away */
  getcurrenttime(&starttime,argv[0]);

  /* Try to run this process with realtime priority */
  if (REALTIME) {
    if(set_realtime_priority() != 0) {
      fprintf(stderr, "%s: Couldn't set realtime priority\n",argv[0]);
      fprintf(stderr, "(You probably need to run this process as root)\n");
    } else
      if (mlockall(MCL_CURRENT|MCL_FUTURE) < 0) {
	fprintf(stderr,"\nError locking memory for realtime\n");
	perror(argv[0]);
      }
  }

  /* Open the comedi device and configure DIO lines for output */
  if (COMEDI) {
    dev = comedi_open(comdevicename);
    if (dev == NULL) {
      comedi_perror(comdevicename);
      fprintf(stderr,"%s: error opening device %s\n",argv[0],comdevicename);
      exit(1);
    }
    /* Check to make sure the subdevice is of the right type */
    stype = comedi_get_subdevice_type(dev,subdevice);
    if (stype != COMEDI_SUBD_DIO) {
      fprintf(stderr,"%s: %d is not a digital I/O subdevice\n",argv[0],subdevice);
      exit(1);
    }
    writemask = 0;
    for (i = 0; i < ndiochannels; i++) {
      writemask = (writemask | (1<<i));
      if (comedi_dio_config(dev,subdevice,i,COMEDI_OUTPUT) < 0) {
	fprintf(stderr,"%s: error configuring digital channels for output\n",argv[0]);
	exit(1);
      }
    }
  }

  /* Read the formatted input, all at once.  That way we know whether
     there are any errors in the input ahead of time. */
  scanret = scanf("%d",&ntransitions);
  if (scanret == 0) {
    fprintf(stderr,"%s: input format not correct\n",argv[0]);
    exit(1);
  } else if (scanret == EOF) {
    fprintf(stderr,"%s: no input!\n",argv[0]);
    exit(1);
  }
  valveseq = (int *) malloc(ntransitions*sizeof(int));
  timeseq = (float *) malloc(ntransitions*sizeof(int));
  if (valveseq == NULL || timeseq == NULL) {
    fprintf(stderr,"%s: error allocating space to parse input\n");
    exit(1);
  }
  for (i = 0; i < ntransitions; i++)
    if ((scanret = scanf("%d %f",valveseq+i,timeseq+i)) == 0 ||
	scanret == EOF) {
      fprintf(stderr,"%s: error parsing input on transition %d\n",argv[0],i);
      exit(1);
    }

  /* Output the start time */
  printf("%dl\t%dl\n",starttime.tv_sec,starttime.tv_usec);
  fflush(stdout);

  /* OK, now we can start processing the valve sequence */
  for (i = 0; i < ntransitions; i++) {
    valvenum = valveseq[i];
    timeswitch = timeseq[i];
    getcurrenttime(&curtime,argv[0]);
    ctimeoffset = (curtime.tv_sec - starttime.tv_sec) +
      (double) (curtime.tv_usec - starttime.tv_usec)/1e6;
    sleeptime = timeswitch - ctimeoffset;
    /* Do we need to wait before changing the valve? */
    if (sleeptime > 0) {
      usleep((unsigned long) (1e6 * sleeptime));
    }
    /* Now manipulate the valve. Since we aren't guaranteed that
       the process sleeps for exactly the time we wanted, get
       the current time again for output */
    getcurrenttime(&curtime,argv[0]);
    ctimeoffset = (curtime.tv_sec - starttime.tv_sec) +
      (double) (curtime.tv_usec - starttime.tv_usec)/1e6;
    if (COMEDI) {
      savevalvenum = valvenum;  /* Prevent read contamination */
      if (comedi_dio_bitfield(dev,subdevice,writemask,&savevalvenum) < 0) {
	fprintf(stderr,"%s: error writing the bitfield to comedi device\n",argv[0]);
	exit(1);
      }
    }
    /* Write output: valve time */
    printf("%d\t%4.4f\n",valvenum,ctimeoffset);
    fflush(stdout);   /* Don't buffer the write, output now */
  }

  if (COMEDI)
    if (comedi_close(dev) < 0) {
      fprintf(stderr,"%s: error closing comedi device\n",argv[0]);
      exit(1);
    }

  free(timeseq);
  free(valveseq);
  exit(0);  /* Everything went according to plan */
}

void getcurrenttime(struct timeval *tvP,char *errmsg)
{
  struct timezone tz;

  if (gettimeofday(tvP,&tz) < 0) {
    perror(errmsg);
    exit(1);
  }
}



int set_realtime_priority(void)
{
  struct sched_param schp;
  /*
   * set the process to realtime privs
   */
  memset(&schp, 0, sizeof(schp));
  schp.sched_priority = sched_get_priority_max(SCHED_FIFO);
  
  if (sched_setscheduler(0, SCHED_FIFO, &schp) != 0) {
    perror("sched_setscheduler");
    return -1;
  }
  
  return 0;
}
