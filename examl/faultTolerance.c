#include <string.h>


/* #include "axml.h" */
#include "globalVariables.h"
#include "faultTolerance.h"
#include "tree.h"



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


/** 
    Obtains a repaired MPI communicator, an updated mapping and a list
    of failed procecsses.    

    @param const tree *tr :                  the main tree structure, will not be modified at all 
    @param List **failedProcs:               an integer list of failed processes
    @param examl_MPI_State **newState:       an intermediate mpi_state,  that will replace our main communicator later, if everything worked out
    @param int **newMapping :                a potential replacement for the mapping. Open work is not assign yet! 

    @return an MPI error code. 

*/
static int getRepairedMPI_State(const tree *tr, examl_MPI_State **newState)
{  
  int
    mpiErr = MPI_SUCCESS ; 

  assert(MPI_SUCCESS == 0); 
  assert(MPI_UNDEFINED < 0 ) ;  /* also just an assumption for convenience */
  
  *newState = calloc(1,sizeof(examl_MPI_State));

  mpiErr |= OMPI_Comm_revoke(mpiState.comm);
  mpiErr |= OMPI_Comm_shrink(mpiState.comm, &((*newState)->comm)); 

  mpiErr |= MPI_Comm_set_errhandler((*newState)->comm, MPI_ERRORS_RETURN); 
  mpiErr |= MPI_Comm_rank((*newState)->comm, &((*newState)->rank)); 
  mpiErr |= MPI_Comm_size((*newState)->comm, &((*newState)->commSize)); 

  return mpiErr;   
}
 


static int getFailedProcesses(int *numFailedRanks, int **failedRanks) 
{
  int 
    i,
    mpiErr = MPI_SUCCESS, 
    failCtr = 0,
    *procArrayPtr = calloc(mpiState.commSize, sizeof(int)),
    *procsTranslatedPtr = calloc(mpiState.commSize, sizeof(int)),
    numFail = 0; 

  MPI_Group failures ; 
  MPI_Group allProc ; 

  for(i = 0; i < mpiState.commSize; ++i)
    procArrayPtr[i] = i; 

  mpiErr |= OMPI_Comm_failure_ack(mpiState.comm); 
  mpiErr |= OMPI_Comm_failure_get_acked(mpiState.comm, &failures);
  mpiErr |= MPI_Comm_group(mpiState.comm, &allProc); 
  mpiErr |= MPI_Group_translate_ranks(allProc, mpiState.commSize, procArrayPtr, failures, procsTranslatedPtr);  
  mpiErr |= MPI_Group_size(failures, &numFail);   

  *failedRanks = calloc(numFail, sizeof(int));   
  int* rkPtr = *failedRanks; 
  for(i = 0; i < mpiState.commSize; ++i)
    if(procsTranslatedPtr[i] != MPI_UNDEFINED) 
      rkPtr[failCtr++] = i ; 

  assert(numFail > 0); 
  assert(numFail == failCtr);

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

  free(newState); 
}



int getMapOldToNewRanks(int *length, int **newRanks, examl_MPI_State *newState)
{
  int
    i,
    *oldRanks = (int*) calloc(newState->commSize, sizeof(int)) ,
    mpiErr = MPI_SUCCESS;
  
  *newRanks = (int*) calloc(mpiState.commSize,sizeof(int)); 
  memset(*newRanks, MPI_UNDEFINED,mpiState.commSize* sizeof(int)); 
  *length = mpiState.commSize; 

  mpiErr |= MPI_Allgather(&(mpiState.rank), 1, MPI_INT, oldRanks, 1, MPI_INT, newState->comm);
  
  /* invert the array  */
  for(i = 0; i < newState->commSize ; ++i) 
    (*newRanks)[oldRanks[i]] = i; 
  
  free(oldRanks);     
  return mpiErr;
}




static void remapWorkAssignment(const tree* const tr, int *length, int **newAssignment, examl_MPI_State *newState)
{
  int
    i,    
    mapLength  = 0,
    *mappedRanks = NULL ;   

  MPI_Group allProc,
    allProcNew ; 
  
  MPI_Comm_group(mpiState.comm, &allProc);
  MPI_Comm_group(newState->comm, &allProcNew);

  getMapOldToNewRanks(&mapLength, &mappedRanks, newState);

  if(tr->manyPartitions)
    {
      pInfo** unassignedPartitions; 
      int 
	ctr = 0, 
	*sitesPerProcess = calloc(newState->commSize,sizeof(int)), 
	unassigned  = 0;
      
      *length = tr->NumberOfModels; 
      *newAssignment = calloc(*length, sizeof(int));

      /* update assignment with new ranks */
      for(i = 0; i < tr->NumberOfModels; ++i ) 
	{
	  int rankNow = mappedRanks[tr->partitionAssignment[i]]; 
	  (*newAssignment)[i] = rankNow; 
	  if(rankNow == MPI_UNDEFINED)
	    ++unassigned; 
	  else 
	    sitesPerProcess[rankNow] += tr->partitionData[i].width ; 	  
	}
      
      unassignedPartitions = malloc(unassigned * sizeof(pInfo*)); 
      for(i = 0; i < tr->NumberOfModels; ++i)
	if((*newAssignment)[i] == MPI_UNDEFINED)
	  unassignedPartitions[ctr++] = &(tr->partitionData[i]);

      qsort(unassignedPartitions,  unassigned, sizeof(pInfo*), partComparePtr); 
      
      /* TODO continue */

      
      free(unassignedPartitions);
      free(sitesPerProcess); 
    }
  else 
    {
      assert(0); 
    }  


  free(mappedRanks); 
}







static int computeMissingSamePhaseSameGen(const tree *tr, tree **localTree, examl_MPI_State *newState)
{
  int 
    length = 0, 
    *newPartitionAssignment = NULL, 
    mpiErr = MPI_SUCCESS;   

  *localTree = calloc(1,sizeof(tree));   
  
  /* setupTree(*localTree);  */

  
  remapWorkAssignment(tr,&length, &newPartitionAssignment, newState); 

  /* TODO  */

  return mpiErr; 
}


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
    *phaseOfProcPtr = NULL,
    *genOfProcPtr = NULL; 
  
  boolean
    erroneousState = TRUE; 
  tree
    *localTree = NULL; 

  examl_MPI_State
    *newState = NULL; 


  if( NOT tr->manyPartitions)
    {
      puts("something went wrong, but FT not implemented for cyclic distribution."); 
      MPI_Abort(mpiState.comm, -1); 
    }    

  puts("something went wrong. Starting recovery.") ; 

  while(erroneousState)
    {
      boolean samePhase = TRUE, 
	sameGen = TRUE; 
      
      erroneousState = FALSE; 

      erroneousState |= getRepairedMPI_State(tr, &newState); 

      int *failedProcs = NULL; 
      int numFailed = 0; 
      erroneousState |= getFailedProcesses(&numFailed,&failedProcs);

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
	  tree *localTree = NULL; 
	  computeMissingSamePhaseSameGen(tr, &localTree, newState);	  
	}
      else 
	{
	  /* TODO  */
	  puts("FT not implemented yet for failure during communication (where some processes got the correct result)."); 
	  assert(0);
	}



      /* TODO repeat problematic communication */
      
      consolidateData(tr, localTree);
      freeLocalTree(localTree); 

      /* free list, if we have to redo it */
      free(failedProcs); 
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
