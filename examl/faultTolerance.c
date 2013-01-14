#include "axml.h"
#include "globalVariables.h"
#include "faultTolerance.h"


/* TODO debug the SEGFAULT */

/* 
   IMPORTANT assumption (not fully verified): ompi_shrink maps each
   rank r to a rank r' = r-f, where f is the number of failed ranks
   with a rank < r
*/

/* 
   Assumption: MPI_SUCCESS == 0 
   just assuming this for development speed; 
   Please replace with proper error handling, if needed. 
 */



#ifdef _USE_RTS

static void appendToFailedProcList(int rank, List** failedProcs)
{
  int *tmp = calloc(1,sizeof(int)); 
  List *elem = calloc(1,sizeof(List)); 
  *tmp = rank; 
  elem->value = tmp ;     
  elem->next = *failedProcs; 
  *failedProcs = elem; 
}



static void freeList(List *elem)
{ 
  List *iter = elem; 
  while(iter != NULL)
    {
      List *tmp = iter->next; 
      free(iter); 
      iter = tmp; 
    } 
}

/** 
    Obtains a repaired MPI communicator, an updated mapping and a list
    of failed procecsses.    

    @param const tree *tr :                  the main tree structure, will not be modified at all 
    @param List **failedProcs:               an integer list of failed processes
    @param examl_MPI_State **newState:       an intermediate mpi_state,  that will replace our main communicator later, if everything worked out
    @param int **newMapping :                a potential replacement for the mapping. Open work is not assign yet! 

    @return an MPI error code. 

*/
static int getRepairedMPI_State(const tree *tr, List **failedProcs, examl_MPI_State **newState, int **newMapping)
{
  int
    mpiErr = MPI_SUCCESS, 
    i, 
    numFail = 0,
    failCtr = 0,
    *procArrayPtr = calloc(mpiState.commSize, sizeof(int)),
    *procsTranslatedPtr = calloc(mpiState.commSize, sizeof(int)); 
  MPI_Group failures ; 
  MPI_Group allProc ; 

  assert(MPI_SUCCESS == 0); 
  assert(MPI_UNDEFINED < 0 ) ;  /* also just an assumption for convenience */
  
  *newState = calloc(1,sizeof(examl_MPI_State));
  *newMapping = calloc(tr->manyPartitions ? tr->NumberOfModels : tr->originalCrunchedLength ,sizeof(int)); 

  for(i = 0; i < mpiState.commSize; ++i)
    procArrayPtr[i] = i; 
  
  mpiErr |= OMPI_Comm_failure_ack(mpiState.comm); 
  mpiErr |= OMPI_Comm_failure_get_acked(mpiState.comm, &failures);
  mpiErr |= MPI_Comm_group(mpiState.comm, &allProc);
  mpiErr |= MPI_Group_translate_ranks(allProc, mpiState.commSize, procArrayPtr, failures, procsTranslatedPtr);   

  for(i = 0; i < mpiState.commSize; ++i)
    {
      if(procsTranslatedPtr[i] != MPI_UNDEFINED)	
	{
	  appendToFailedProcList(i, failedProcs);
	  failCtr++;
	}
    }

  /* map sites to new rank system */
  {
    int
      *pos = tr->manyPartitions ? tr->partitionAssignment : tr->siteToProcess,
      *posEnd = tr->manyPartitions ? tr->partitionAssignment + tr->NumberOfModels : tr->siteToProcess + tr->originalCrunchedLength; 
    int ctr = 0; 
    for(; pos != posEnd; ++pos)
      {
	if(procsTranslatedPtr[*pos] == MPI_UNDEFINED)
	  printf("[%d] %d -> OPEN\n", mpiState.rank, *pos); 
	else 
	  printf("[%d] %d -> %d\n",mpiState.rank, *pos, procsTranslatedPtr[*pos] ); 
	(*newMapping)[ctr++] = procsTranslatedPtr[*pos]; /* NOTE: open work is assign MPI_UNDEFINED */
      }
  }
  
  mpiErr |= MPI_Group_size(failures, &numFail); 
  
  assert(numFail != 0); 
  assert(numFail == failCtr); 
   
  mpiErr |= OMPI_Comm_revoke(mpiState.comm);
  mpiErr |= OMPI_Comm_shrink(mpiState.comm, &((*newState)->comm)); 

  mpiErr |= MPI_Comm_set_errhandler((*newState)->comm, MPI_ERRORS_RETURN); 
  mpiErr |= MPI_Comm_rank((*newState)->comm, &((*newState)->rank)); 
  mpiErr |= MPI_Comm_size((*newState)->comm, &((*newState)->commSize)); 
  mpiErr |= MPI_Group_free(&failures); 
  
  free(procArrayPtr); 
  free(procsTranslatedPtr); 

  return mpiErr;   
}
 

static void restoreComm(examl_MPI_State *newState)
{  
  mpiState.comm = newState->comm; 
  mpiState.rank = newState->rank; 
  mpiState.commSize = newState->commSize; 
  mpiState.mpiError = MPI_SUCCESS; 

  /* :TODO: is that all? */

  free(newState); 
}


/* void initLocalTreeForFailedProcs(tree *tr, tree *laocalTree, List *failedProcs) */
/* { */
/*   assert(0);  */
/* } */



void execForFailedProcs(tree *tr)
{  
  assert(0); 
}



void freeLocalTree(tree *localTree)
{
  assert(0);
}



void consolidateData(tree *tr, tree *localTree)
{
  assert(0); 


  
  freeLocalTree(localTree); 
}

/** 
    Top-level function for handling MPI errors. 

    Tries to achieve a state that equals the state from which this
    function was called, but where no error occurred.
    
    A side effect of this function is that sites/partitions may be re-distributed. 
    
    @param tree *tr:                 the main tree structure 
*/
void handleMPIError(tree *tr)
{  
  int
    i,
    *newMapping = NULL,  	/* of either the sites or the partitions to processes */
    *phaseOfProcPtr = NULL,
    *genOfProcPtr = NULL; 
  
  boolean
    erroneousState = TRUE; 
  tree
    *localTree = NULL; 

  examl_MPI_State
    *newState = NULL; 

  while(erroneousState)
    {
      boolean samePhase = TRUE, 
	sameGen = TRUE; 
      
      erroneousState = FALSE; 

      List *failedProcsUnprocessed = NULL ;      
      erroneousState |= getRepairedMPI_State(tr, &failedProcsUnprocessed, &newState, &newMapping); 
      
      /* notify */
      List *iter; 
      for(iter = failedProcsUnprocessed; iter != NULL; iter = iter->next)
	printf("%d detected failure of %d\n", mpiState.rank, ACCESS_LIST_INT(iter)); 
      
      phaseOfProcPtr = calloc(newState->commSize,  sizeof(int)); 
      genOfProcPtr = calloc(newState->commSize, sizeof(int)); 
      
      /* determine phase and generation of other processes  */ 
      erroneousState |= MPI_Allgather(&(mpiState.commPhase), 1, MPI_INT, phaseOfProcPtr,1,MPI_INT,newState->comm);
      erroneousState |= MPI_Allgather(&(mpiState.generation[mpiState.commPhase]), 1, MPI_INT, genOfProcPtr, 1, MPI_INT, newState->comm); 
      
      int onePhase = phaseOfProcPtr[0],
	otherPhase = MPI_UNDEFINED, 
	oneGen = genOfProcPtr[0],
	otherGen = MPI_UNDEFINED; 
      for(i = 0; i < newState->commSize; ++i)
	{
	  samePhase &= (onePhase == phaseOfProcPtr[i]);
	  if(NOT samePhase)
	    {
	      if(otherPhase == MPI_UNDEFINED)
		otherPhase = phaseOfProcPtr[i]; 
	      else if(otherPhase != phaseOfProcPtr[i] && onePhase != phaseOfProcPtr[i])
		{
		  printf("something seriously went wrong during fault recovery:  more than 3 phases (%d,%d,%d).\nGiving up, this is a bug.\n", onePhase, otherPhase, phaseOfProcPtr[i]); 
		  assert(0); 
		}
	    }


	  sameGen &= (oneGen == genOfProcPtr[i]); 
	  if(NOT sameGen)
	    {
	      if(otherGen == MPI_UNDEFINED)
		otherGen = genOfProcPtr[i]; 
	      else if(otherGen != genOfProcPtr[i] && oneGen != genOfProcPtr[i])
		{
		  printf("something seriously went wrong during fault recovery:  more than 3 generations (%d,%d,%d).\nGiving up, this is a bug.\n", oneGen, otherGen, genOfProcPtr[i]); 
		  assert(0); 
		}
	    }	  
	}
      
      if(sameGen && samePhase)
	{ 
	  /* TODO  */
	  puts("not implemented yet."); 
	  assert(0); 
	}
      else 
	{
	  /* TODO  */
	  puts("not implemented yet."); 
	  assert(0);
	}

      /* free list, if we have to redo it */
      freeList(failedProcsUnprocessed);
      free(phaseOfProcPtr); 
      free(genOfProcPtr);
      freeLocalTree(localTree); 
    }
  
  /* if everything worked out  */
  restoreComm(newState);
  consolidateData(tr, localTree);
}



#else 
void doNothing()
{
}
#endif
