#ifndef _THREAD_H
#define _THREAD_H
 
#include "globalVariables.h"

void startPthreads(int argc, char *argv[]); 

void hybrid_allreduce(tree *tr,  size_t length); 
void hybrid_allreduce_evaluate(tree *tr, size_t length);
void hybrid_allreduce_makenewz(tree *tr, size_t length);

void tb_threadsWait(tree *tr);
void tb_lockThreads(tree *tr) ;
void tb_unlockThreads(tree *tr); 
void tb_workerTrap(tree *tr);



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
	tb_unlockThreads(tr);						\
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
	tb_unlockThreads(tr);						\
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
	tb_unlockThreads(tr);					\
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
