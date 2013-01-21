#ifdef _HYBRID

/* :TODO: we could migrate to the c11 standard libraries threads at
   some point. It may still take a while until this is
   implemented.  */
/* #ifdef __STDC_NO_THREADS__ */
#include <pthread.h>
/* #else  */
/* #include <threads.h> */
/* #endif */


#include "axml.h"
#include "globalVariables.h"

/* 
   :TODO: stuff I omitted for simplicity 

   * Thread to core pinning 

*/

/* extern int realMain(tree *tr, analdef *adef, int tid) ; */
extern int realMain(int tid, int argc, char *argv[]); 


static void* realMainThreadWrapper(void *tmp)
{
  threadData *tData = (threadData*) tmp; 

  realMain(tData->tid ,tData->argc, tData->argv); 

  pthread_exit(NULL); 
  return NULL; 
}


/**
   @brief Starts the pthreads. Modifies the mpi singleton accordingly.
 */
void startPthreads(int argc, char *argv[])
{
  pthread_attr_t attr;
  int rc, t;
  threadData *tData;

  /* jobCycle        = 0; */
  /* threadJob       = 0; */

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

  /* pthread_mutex_init(&mutex , (pthread_mutexattr_t *)NULL); */

  mpiState.threads    = (pthread_t *)malloc(mpiState.numberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)malloc(mpiState.numberOfThreads * sizeof(threadData));
  /* reductionBuffer          = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels); */
  /* reductionBufferTwo       = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels); */
  /* reductionBufferThree     = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels); */
  /* reductionBufferParsimony = (volatile int *)malloc(sizeof(volatile int) *  NumberOfThreads); */

  /* mpiState.barrier = malloc(sizeof(volatile char) *  mpiState.numberOfThreads); */

  /* for(t = 0; t < mpiState.numberOfThreads; t++) */
    /* mpiState.barrier[t] = 0; */

  /* branchInfos              = (volatile branchInfo **)malloc(sizeof(volatile branchInfo *) * NumberOfThreads); */

  for(t = 1; t < mpiState.numberOfThreads; t++)
  {
    tData[t].tid  = t;
    tData[t].argc  = argc;
    tData[t].argv  = argv;

    rc = pthread_create(&(mpiState.threads[t]), &attr, realMainThreadWrapper, (void *)(&tData[t]));
    if(rc)
      {
	printf("ERROR; return code from pthread_create() is %d\n", rc);
	exit(-1);
      }
  }
}



void threadBarrier(int tid)
{
  int i;

  if(tid == 0)
    {
      for(i = 0 ;i < mpiState.numberOfThreads; ++i)
        mpiState.barrier[i] = FALSE;
      mpiState.barrierIsCrossed = FALSE;
      mpiState.threadsCanCheckBarrier = TRUE;
    }
  else
    while(NOT mpiState.threadsCanCheckBarrier);

  mpiState.barrier[tid] = TRUE;

  if(tid == 0)
    {
      nat sum;
      do
        {
          sum = 0;
          for(i = 0; i < mpiState.numberOfThreads; ++i)
            sum += mpiState.barrier[i];
        }
      while(sum < mpiState.numberOfThreads);

      mpiState.threadsCanCheckBarrier = FALSE;
      mpiState.barrierIsCrossed = TRUE;
    }
  else
    while(NOT mpiState.barrierIsCrossed);
}


#else 
#include "axml.h"
void dummyFunction(int a); 

void dummyFunction(int a)
{
  printf("o my god, this is annoying %d.\n" ,a ); 
}


void threadBarrier(int tid)
{
}
#endif
