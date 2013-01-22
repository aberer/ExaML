


void threadBarrier(int tid); 

#ifdef _HYBRID
#define ABS_ID(tr) ((mpiState.rank  * mpiState.numberOfThreads) + tr->threadId)
#define ABS_NUM_RANK (mpiState.commSize *  mpiState.numberOfThreads)


#define MASTER_TREE (mpiState.allTrees[0] )
#define GET_TREE_NUM(x) (mpiState.allTrees[x])


/* TODO can i trust this barrier?  */
/* :TODO: not very efficient, is it?  */

/** @notice only works with sums, but that is okay  */
#define HYBRID_ALLREDUCE_VAR(tr, tree_var, length, mpi_type, type)	\
  {									\
    if(tr->threadId == 0)						\
      {									\
	int i,j;							\
	type *buf = calloc(length, sizeof(type)) ;			\
									\
	threadBarrier(tr->threadId);					\
  									\
	for(i = 0; i < length;++i)					\
	  for(j = 1 ; j < mpiState.numberOfThreads ; ++j)		\
	    tr->tree_var[i] +=  GET_TREE_NUM(j)->tree_var[i]; \
  									\
	MPI_Allreduce(tr->tree_var, buf, length, mpi_type, MPI_SUM, mpiState.comm); \
  									\
	for(j = 0; j < mpiState.numberOfThreads; ++j)			\
	  memcpy(GET_TREE_NUM(j)->tree_var, buf, sizeof(type) * length); \
									\
	free(buf);							\
	threadBarrier(tr->threadId);					\
      }									\
    else								\
      {									\
	/* ensures that the threads have provided their results  */	\
	threadBarrier(tr->threadId);					\
      									\
	/* ensures that the threads received their results */		\
	threadBarrier(tr->threadId);					\
      }									\
  }									\


/* TODO 
 * with error later
 */
#define HYBRID_BCAST_VAR(tr,tree_var,length,mpi_type,type)		\
  {									\
    if(tr->threadId == 0)						\
      {									\
	MPI_Bcast(tr->tree_var,length,mpi_type,0,mpiState.comm);	\
	threadBarrier(tr->threadId);					\
      }									\
    else								\
      {									\
	threadBarrier(tr->threadId);					\
	memcpy(tr->tree_var, MASTER_TREE->tree_var , length * sizeof(type) ); \
      }									\
  }									\



#define HYBRID_BCAST_VAR_1(tr,tree_var,mpi_type)			\
  {									\
   if(tr->threadId ==0)							\
      {									\
	MPI_Bcast(&(tr->tree_var),1,mpi_type,0,mpiState.comm);		\
	threadBarrier(tr->threadId);					\
      }									\
    else								\
      {									\
	threadBarrier(tr->threadId);					\
	tr->tree_var = MASTER_TREE->tree_var;				\
      }									\
  }									\


#define HYBRID_BARRIER(tid)			\
  {						\
   threadBarrier(tid);				\
   if(tid == 0)					\
     MPI_Barrier(mpiState.comm);		\
  }						\


#else 

#define ABS_ID(tr) (mpiState.rank)
#define ABS_NUM_RANK (mpiState.commSize)


#define HYBRID_BCAST_VAR(tr, tree_var,length, mpi_type,type)	\
  {								\
    MPI_Bcast(tr->tree_var, length, mpi_type,0,mpiState.comm);	\
  }								\

#define HYBRID_BCAST_VAR_1(tr,tree_var,mpi_type)	\
  {							\
    MPI_Bcast(&(tr->tree_var), 1,mpi_type,0,mpiState.comm);	\
  }							\

#define HYBRID_BARRIER(tid)			\
  {MPI_Barrier(mpiState.comm);}			\


#define HYBRID_ALLREDUCE_VAR(tr, tree_var, length, mpi_type, type)	\
  {									\
    type *buf = calloc(length, sizeof(type));				\
    MPI_Allreduce(tr->tree_var, buf, length, mpi_type, MPI_SUM, mpiState.comm); \
    memcpy(tr->tree_var, buf, length * sizeof(type));			\
    free(buf);								\
  }									\

#endif
