#ifdef _USE_PTHREADS

#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>

#include "axml.h"
#include "globalVariables.h"


/* 
   :TODO: stuff I omitted for simplicity 

   * Thread to core pinning 



*/



/* static void execFunction(tree *tr, tree *localTree, int tid, int n) */
/* { */
/*   int */
/*     currentJob */
/*     ; */

/*   currentJob = threadJob >> 16; */

/*   switch(currentJob) */
/*   {             */
/*   case 0:  */
/*   case 13:  */
/*     { */
/*       printf("worker %d/%d has nothing to do\n", tid, processID);  */
/*       sleep(1); */
/*       break;  */
/*     } */
/*   default: */
/*     printf("Job %d\n", currentJob); */
/*     assert(0); */
/*   } */

  
/* } */



/* void masterBarrier(int jobType, tree *tr) */
/* {   */
/*   const int  */
/*     n = NumberOfThreads; */

/*   int  */
/*     i,  */
/*     sum; */

/*   jobCycle = !jobCycle; */
/*   threadJob = (jobType << 16) + jobCycle; */

/*   /\* execFunction(tr, tr, 0, n); *\/ */


/*   do */
/*   { */
/*     for(i = 1, sum = 1; i < n; i++) */
/*       sum += barrierBuffer[i]; */
/*   } */
/*   while(sum < n);   */

/*   for(i = 1; i < n; i++) */
/*     barrierBuffer[i] = 0; */
/* } */ 



/* static void *likelihoodThread(void *tData) */
/* { */
/*   threadData *td = (threadData*)tData; */
/*   tree */
/*     *tr = td->tr, */
/*     *localTree = (tree *)malloc(sizeof(tree)); */
/*   int */
/*     myCycle = 0; */

/*   const int  */
/*     n = NumberOfThreads, */
/*       tid = td->threadNumber; */

/* /\* #ifndef _PORTABLE_PTHREADS *\/ */
/*   /\* :TODO: *\/ */
/*   /\* pinToCore(tid); *\/ */
/* /\* #endif *\/ */

/*   printf("\n\tThis is the ExaML Worker Pthread: %d\tprocess rank: %d\n", tid, processID); */

/*   while(1) */
/*   { */
/*     while (myCycle == threadJob); */
/*     myCycle = threadJob; */

/*     execFunction(tr, localTree, tid, n); */


/*     barrierBuffer[tid] = 1;      */
/*   } */

/*   return (void*)NULL; */
/* } */


/* void startPthreads(tree *tr) */
/* { */
/*   pthread_t *threads; */
/*   pthread_attr_t attr; */
/*   int rc, t; */
/*   threadData *tData; */

/*   jobCycle        = 0; */
/*   threadJob       = 0; */

/*   printf("\n\tThis is the ExaML Master Pthread\tprocess rank: %d\n", processID); */

/*   pthread_attr_init(&attr); */
/*   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED); */

/*   pthread_mutex_init(&mutex , (pthread_mutexattr_t *)NULL); */

/*   threads    = (pthread_t *)malloc(NumberOfThreads * sizeof(pthread_t)); */
/*   tData      = (threadData *)malloc(NumberOfThreads * sizeof(threadData)); */

/*   barrierBuffer            = (volatile char *)malloc(sizeof(volatile char) *  NumberOfThreads); */

/*   for(t = 0; t < NumberOfThreads; t++) */
/*     barrierBuffer[t] = 0; */

/*   /\* printf("number of threads is %d\n", NumberOfThreads);  *\/ */


/*   for(t = 1; t < NumberOfThreads; t++) */
/*     { */
/*       tData[t].tr  = tr; */
/*       tData[t].threadNumber = t; */
/*       rc = pthread_create(&threads[t], &attr, likelihoodThread, (void *)(&tData[t])); */
/*       if(rc) */
/* 	{ */
/* 	  printf("ERROR; return code from pthread_create() is %d\n", rc); */
/* 	  exit(-1); */
/* 	} */
/*     } */


/* } */
#else 
#include "axml.h"
void dummyFunction(int a); 

void dummyFunction(int a)
{
  printf("o my god, this is annoying %d.\n" ,a ); 
}

#endif
