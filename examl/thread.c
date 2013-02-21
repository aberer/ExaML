
#include "axml.h"
#include "globalVariables.h"
#include "thread.h"


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




void tb_waitForWorkers(tree *tr) 
{
  assert(threadsAreLocked == FALSE); 
  assert(tr->threadId == 0); 
  /* DM(tr,"master locks workers\n");  */
  mpiState.barrier[0] = 1;  
  int sum;
  do
    {
      sum = 0;
      for(int i = 0; i < mpiState.numberOfThreads; ++i)
	if(mpiState.barrier[i])
	  ++sum; 
    } while(sum != mpiState.numberOfThreads);
  threadsAreLocked = TRUE;   
}


void tb_workerWait(tree *tr)
{
  assert(tr->threadId > 0); 
  /* DM(tr, " is blocked\n") ; */
  mpiState.barrier[tr->threadId] = TRUE; 
  while(mpiState.barrier[tr->threadId]); 
}


void tb_unlockThreads(tree *tr)
{
  assert(threadsAreLocked == TRUE); 
  assert(tr->threadId == 0);
  memset((void*)mpiState.barrier, 0, sizeof(volatile char) * mpiState.numberOfThreads);
  /* for(int i = 0; i < mpiState.numberOfThreads; ++i) */
  /*   mpiState.barrier[i] = FALSE; */
  threadsAreLocked = FALSE; 


  /* DM(tr, "master unlocked threads\n");  */
}


void tb_barrier(tree *tr)
{
  if(tr->threadId == 0 ) 
    tb_waitForWorkers(tr); 
  else 
    tb_workerWait(tr); 
}


void hybrid_allreduce_makenewz(tree *tr, size_t length)
{
  if(tr->threadId == 0)
    {
      tb_waitForWorkers(tr); 

      for(int t = 1; t < mpiState.numberOfThreads; ++t)
	{
	  double *otherPtr = GET_TREE_NUM(t)->reductionBuffer; 
	  
	  for(int j = 0; j < length; ++j)
	    tr->reductionBuffer[j] += otherPtr[j]; 
	}
	
      MPI_Allreduce(MPI_IN_PLACE, tr->reductionBuffer, length, MPI_DOUBLE, MPI_SUM, mpiState.comm); 
      tb_unlockThreads(tr); 
      tb_waitForWorkers(tr); 
      tb_unlockThreads(tr); 
    }
  else 
    {
      tb_workerWait(tr); 
      memcpy(tr->reductionBuffer, MASTER_TREE->reductionBuffer, sizeof(double) * length);
      tb_workerWait(tr); 
    }
}



void hybrid_allreduce_evaluate(tree *tr, size_t length)
{
  /* only the local reduction  */
  if(numRounds != 0)
    {
      int rounds = 0; 
      int dist = 1;
      
      while(rounds < numRounds)
	{
	  tb_barrier(tr); 
	  if(tr->threadId == 0)
	    tb_unlockThreads(tr); 

	  int otherId = tr->threadId + dist; 
	  if( ( (  tr->threadId & ((dist << 1) - 1 )  ) == 0 ) 
	      &&  (otherId  < mpiState.numberOfThreads) )
	    {
	      double *otherPtr = GET_TREE_NUM(otherId)->perPartitionLH; 

	      for(int i = 0; i < length ;++i)
		tr->perPartitionLH[i] += otherPtr[i];
	    }

	  dist <<= 1 ; 
	  ++rounds; 
	}
    }

  /* mpi reduce and local broadcast */
  if(tr->threadId == 0)
    {
      tb_waitForWorkers(tr);
      MPI_Allreduce(MPI_IN_PLACE, tr->perPartitionLH, length, MPI_DOUBLE, MPI_SUM, mpiState.comm); 
      tb_unlockThreads(tr); 

      tb_waitForWorkers(tr); 
      tb_unlockThreads(tr); 
    }
  else 
    {
      tb_workerWait(tr); 
      memcpy(tr->perPartitionLH, MASTER_TREE->perPartitionLH, sizeof(double) * length);
      tb_workerWait(tr); 
    }
}




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

  mpiState.barrier = (volatile char*) calloc(mpiState.numberOfThreads, sizeof(volatile char)); 
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
  }
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
