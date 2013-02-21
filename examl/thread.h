#ifndef _THREAD_H
#define _THREAD_H
 
#include "globalVariables.h"

void startPthreads(int argc, char *argv[]); 

void hybrid_allreduce(tree *tr,  size_t length); 
void hybrid_allreduce_evaluate(tree *tr, size_t length);
void hybrid_allreduce_makenewz(tree *tr, size_t length);

void tb_workerTrap_worker(tree *tr); 
void tb_workerTrap_master(tree *tr) ; 



inline void tb_releaseWorkers(tree *tr)
{  
  ++mpiState.globalGen; 
  mpiState.threadsAreLocked = FALSE ; 
}


inline void tb_workerTrap(tree *tr)
{
  if(tr->threadId == 0 ) 
    tb_workerTrap_master(tr); 
  else 
    tb_workerTrap_worker(tr); 
}


/* :TODO: god this inlining sucks, figure this out later  */
#ifdef _DEBUG
int ABS_ID(int num); 
int ABS_NUM_RANK(); 
#else 
inline int ABS_ID(int num)
{
  assert(mpiState.numberOfThreads > 0); 
  return(mpiState.rank * mpiState.numberOfThreads) + num ; 
}


inline int ABS_NUM_RANK()
{
  assert(mpiState.numberOfThreads > 0); 
  return mpiState.commSize * mpiState.numberOfThreads; 
}
#endif

#ifdef _HYBRID


#define MASTER_TREE (mpiState.allTrees[0] )
#define GET_TREE_NUM(x) (mpiState.allTrees[x])


#define HYBRID_ALLREDUCE_VAR_TWO_BARRIERS_TREE(tr, tree_var, length, mpi_type, type) \
  {									\
  int rounds = 0;							\
  int dist = 1;								\
  while(rounds < numRounds)						\
    {									\
									\
      int potRecv = ( tr->threadId - dist );				\
      if( potRecv %  (dist << 1 )  == 0 )				\
	{								\
	  /* DM(tr, " sending to %d in round %d, distance was %d\n", potRecv,rounds, dist );  */ \
	  for(int i = 0; i < length; ++i)				\
	    GET_TREE_NUM(potRecv)->tree_var[i] += tr->tree_var[i];	\
	}								\
									\
      tb_workerTrap(tr);						\
      if(tr->threadId==0)						\
	tb_releaseWorkers(tr);						\
									\
      dist <<= 1 ;							\
      ++rounds;								\
    }									\
									\
  if(tr->threadId == 0)							\
    {									\
      /* DM(tr,"\n"); */							\
      tb_workerTrap(tr);						\
      MPI_Allreduce(MPI_IN_PLACE, tr->tree_var , length, mpi_type, MPI_SUM, mpiState.comm); \
      tb_releaseWorkers(tr);						\
      tb_workerTrap(tr);						\
      tb_releaseWorkers(tr);						\
    }									\
  else									\
    {									\
      tb_workerTrap(tr);						\
      memcpy(tr->tree_var, MASTER_TREE->tree_var, sizeof(type) * length); \
      tb_workerTrap(tr);						\
    }									\
  }									\


/* TODO can i trust this barrier?  */
/* :TODO: not very efficient, is it?  */

/** @notice only works with sums, but that is okay  */
#define HYBRID_ALLREDUCE_VAR_ONE_BARRIER(tr, tree_var, length, mpi_type, type)	\
  {									\
    if(tr->threadId == 0)						\
      {									\
	tb_workerTrap(tr);						\
  									\
	for(int i = 0; i < length;++i)					\
	  for(int j = 1 ; j < mpiState.numberOfThreads ; ++j)		\
	    tr->tree_var[i] +=  GET_TREE_NUM(j)->tree_var[i];		\
					 				\
	MPI_Allreduce(MPI_IN_PLACE, tr->tree_var , length, mpi_type, MPI_SUM, mpiState.comm); \
									\
	for(int i = 1; i < mpiState.numberOfThreads; ++i)		\
	  memcpy(GET_TREE_NUM(i)->tree_var, MASTER_TREE->tree_var,sizeof(type) * length); \
									\
	tb_releaseWorkers(tr);						\
      }									\
    else								\
      tb_workerTrap(tr);						\
  }									\


#define HYBRID_ALLREDUCE_VAR_TWO_BARRIERS(tr, tree_var, length, mpi_type, type)	\
  {									\
    if(tr->threadId == 0)						\
      {									\
	tb_workerTrap(tr);						\
  									\
	for(int i = 0; i < length;++i)					\
	  for(int j = 1 ; j < mpiState.numberOfThreads ; ++j)		\
	    tr->tree_var[i] +=  GET_TREE_NUM(j)->tree_var[i];		\
					 				\
	MPI_Allreduce(MPI_IN_PLACE, (void*)tr->tree_var , length, mpi_type, MPI_SUM, mpiState.comm); \
	tb_releaseWorkers(tr);						\
	tb_workerTrap(tr);						\
	tb_releaseWorkers(tr);						\
      }									\
    else								\
      {									\
	tb_workerTrap(tr);						\
	memcpy((void*)tr->tree_var, (void*)MASTER_TREE->tree_var, sizeof(type) * length); \
	tb_workerTrap(tr);						\
      }									\
  }									\




#define HYBRID_ALLREDUCE_VAR HYBRID_ALLREDUCE_VAR_TWO_BARRIERS_TREE
/* #define HYBRID_ALLREDUCE_VAR HYBRID_ALLREDUCE_VAR_ONE_BARRIER */
/* #define HYBRID_ALLREDUCE_VAR HYBRID_ALLREDUCE_VAR_TWO_BARRIERS */



/* TODO 
 * with error later
 */
#define HYBRID_BCAST_VAR(tr,tree_var,length,mpi_type,type)		\
  {									\
    if(tr->threadId == 0)						\
      {									\
	tb_workerTrap(tr);						\
	MPI_Bcast(tr->tree_var,length,mpi_type,0,mpiState.comm);	\
	for(int i = 1; i < mpiState.numberOfThreads; ++i)		\
	  memcpy(GET_TREE_NUM(i)->tree_var, MASTER_TREE->tree_var , length * sizeof(type) ); \
	tb_releaseWorkers(tr);						\
      }									\
    else								\
      {									\
	tb_workerTrap(tr);						\
      }									\
  }									\


#define HYBRID_BCAST_VAR_1(tr,tree_var,mpi_type)			\
  {									\
   if(tr->threadId ==0)							\
      {									\
	tb_workerTrap(tr);						\
	MPI_Bcast(&(tr->tree_var),1,mpi_type,0,mpiState.comm);		\
	for(int i = 1; i < mpiState.numberOfThreads; ++i)		\
	  GET_TREE_NUM(i)->tree_var = MASTER_TREE->tree_var;		\
	tb_releaseWorkers(tr);						\
      }									\
    else								\
      {									\
	tb_workerTrap(tr);						\
      }									\
  }									\


#define HYBRID_BARRIER(tr)					\
  {								\
    tb_workerTrap(tr);						\
    if(tr->threadId == 0)					\
      {								\
	MPI_Barrier(mpiState.comm);				\
	tb_releaseWorkers(tr);					\
      }								\
  }								\
    

#else 

/* #define ABS_ID(tr) (mpiState.rank) */
/* #define ABS_NUM_RANK (mpiState.commSize) */


#define HYBRID_BCAST_VAR(tr, tree_var,length, mpi_type,type)	\
  {								\
    MPI_Bcast(tr->tree_var, length, mpi_type,0,mpiState.comm);	\
  }								\

#define HYBRID_BCAST_VAR_1(tr,tree_var,mpi_type)	\
  {							\
    MPI_Bcast(&(tr->tree_var), 1,mpi_type,0,mpiState.comm);	\
  }							\

#define HYBRID_BARRIER(tr)			\
  {MPI_Barrier(mpiState.comm);}			\


#define HYBRID_ALLREDUCE_VAR(tr, tree_var, length, mpi_type, type)	\
  {									\
    type *buf = calloc(length, sizeof(type));				\
    MPI_Allreduce(tr->tree_var, buf, length, mpi_type, MPI_SUM, mpiState.comm); \
    memcpy(tr->tree_var, buf, length * sizeof(type));			\
    free(buf);								\
  }									\

#endif


#endif
