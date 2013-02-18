
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
  
  /* printf("\n\n%d enters wrapper\n\n", tData->tid);  */

  realMain(tData->tid ,tData->argc, tData->argv); 

  printf("\n\nAttention: thread %d exits\n\n", tData->tid); 
  pthread_exit(NULL); 
  return NULL; 
}


/**
   @brief Starts the pthreads. Modifies the mpi singleton accordingly.
 */
void startPthreads(int argc, char *argv[])
{
  pthread_attr_t* attr = malloc(sizeof(pthread_attr_t));
  int rc, t;
  threadData *tData;

  pthread_attr_init(attr);
  /* pthread_attr_setdetachstate(attr, PTHREAD_CREATE_DETACHED); */

  mpiState.threads    = (pthread_t *)malloc(mpiState.numberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)malloc(mpiState.numberOfThreads * sizeof(threadData));

  mpiState.exitCodes = calloc(mpiState.numberOfThreads,sizeof(int));

  mpiState.localGen = calloc(mpiState.numberOfThreads, sizeof(int));
  mpiState.globalGen = 0; 
  mpiState.threadsAreLocked = FALSE; 

  mpiState.allTrees = calloc(mpiState.numberOfThreads, sizeof(tree*)); 

  for(t = 1; t < mpiState.numberOfThreads; t++)
  {
    tData[t].tid  = t;
    tData[t].argc  = argc;
    tData[t].argv  = argv;

    rc = pthread_create(&(mpiState.threads[t]), attr, realMainThreadWrapper, (void *)(&tData[t]));

    if(rc)
      {
	printf("ERROR; return code from pthread_create() is %d\n", rc);
	exit(-1);
      } 
    
    while(mpiState.localGen[t] != 1) ; 
  }
}


void tb_workerTrap(tree *tr)
{
  if(tr->threadId == 0)
    {
      assert(NOT mpiState.threadsAreLocked); 
      mpiState.localGen[tr->threadId]++; 
      const int myGen = mpiState.localGen[tr->threadId]; 
      
      int sum;
      do 
	{
	  sum = 0; 
	  for(int i = 0; i < mpiState.numberOfThreads; ++i)
	    if(myGen == mpiState.localGen[i])
	      sum++; 
	} while(sum != mpiState.numberOfThreads); 
      mpiState.threadsAreLocked = TRUE; 
    }
  else 
    {      
      mpiState.localGen[ tr->threadId ]++; 
      const int myGen = mpiState.localGen[tr->threadId]; 
      while(myGen != mpiState.globalGen); 
    }
}



void tb_releaseWorkers(tree *tr)
{  
  assert(tr->threadId == 0 && mpiState.threadsAreLocked); 
  ++mpiState.globalGen; 
  mpiState.threadsAreLocked = FALSE; 
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

#endif
