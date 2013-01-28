
#include "axml.h"
#include "globalVariables.h"


#ifdef _DEBUG
int ABS_ID(int num)
{
  assert(mpiState.numberOfThreads > 0); 
  return  (mpiState.rank * mpiState.numberOfThreads) + num ; 
}


int ABS_NUM_RANK()
{
  assert(mpiState.numberOfThreads > 0); 
  return mpiState.commSize * mpiState.numberOfThreads; 
}

#endif



#ifdef _HYBRID

/* :TODO: we could migrate to the c11 standard libraries threads at
   some point. It may still take a while until this is
   implemented.  */
/* #ifdef __STDC_NO_THREADS__ */
#include <pthread.h>
/* #else  */
/* #include <threads.h> */
/* #endif */






/* 
   :TODO: stuff I omitted for simplicity 

   * Thread to core pinning 
*/

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

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

  mpiState.threads    = (pthread_t *)malloc(mpiState.numberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)malloc(mpiState.numberOfThreads * sizeof(threadData));
  /* pthread_mutex_init(&mutex , (pthread_mutexattr_t *)NULL); */
  /* reductionBuffer          = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels); */
  /* reductionBufferTwo       = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels); */
  /* reductionBufferThree     = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels); */
  /* reductionBufferParsimony = (volatile int *)malloc(sizeof(volatile int) *  NumberOfThreads); */

  /* branchInfos              = (volatile branchInfo **)malloc(sizeof(volatile branchInfo *) * NumberOfThreads); */

  
  mpiState.exitCodes = calloc(mpiState.numberOfThreads,sizeof(int));

  /* mpiState.barrier = calloc(mpiState.numberOfThreads, sizeof(volatile boolean) ); */
  /* for(t = 0; t < mpiState.numberOfThreads; t++) */
  /*   mpiState.barrier[t] = FALSE; */

  mpiState.allTrees = calloc(mpiState.numberOfThreads, sizeof(tree*)); 
  

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
  pthread_barrier_wait(&(mpiState.pBarrier)); 
}


#else 
#include "axml.h"


inline int ABS_ID(int num)
{
  return  (mpiState.rank) + num ; 
}


inline int ABS_NUM_RANK()
{
  return mpiState.commSize; 
}

inline void threadBarrier(int tid)
{
}
#endif
