
#include "axml.h"
#include "globalVariables.h"
#include "thread.h"



int isThisMyPartition(tree *tr, int model)
{  
  assert(tr->manyPartitions); 
  return ( tr->partitionAssignment[model] == ABS_ID(tr) ) ; 
}



static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
    y = 362436069,
    z = 21288629,
    w = 14921776,
    c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}


boolean setupTree (tree *tr)
{
  nodeptr  p0, p, q;
  int
    i,
    j,   
    tips,
    inter; 

  
  tr->bigCutoff = FALSE;
  
  tr->maxCategories = MAX(4, tr->categories);
  
  tr->partitionContributions = (double *)malloc(sizeof(double) * tr->NumberOfModels);
  
  for(i = 0; i < tr->NumberOfModels; i++)
    tr->partitionContributions[i] = -1.0;
  
  tr->perPartitionLH = calloc(tr->NumberOfModels, sizeof(double));

  tips  = tr->mxtips;
  inter = tr->mxtips - 1; 
  
  tr->fracchanges  = (double *)malloc(tr->NumberOfModels * sizeof(double)); 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)calloc(tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)calloc(tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  
            
  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);
  tr->td[0].executeModel = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);
  tr->td[0].parameterValues = (double *)malloc(sizeof(double) * tr->NumberOfModels);
  
  for(i = 0; i < tr->NumberOfModels; i++)
    tr->fracchanges[i] = -1.0;
  tr->fracchange = -1.0;
  
  tr->constraintVector = (int *)malloc((2 * tr->mxtips) * sizeof(int));
  
  tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));
   

  if (!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }
  
  tr->nodeBaseAddress = p0;


  if (!(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
    {
      p = p0++;

      p->hash   =  KISS32(); /* hast table stuff */
      p->x      =  0;
      p->xBips  =  0;
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL;     
      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++)
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++)
	{	 
	  p = p0++;
	  if(j == 1)
	    {
	      p->xBips = 1;
	      p->x = 1;
	    }
	  else
	    {
	      p->xBips = 0;
	      p->x =  0;
	    }
	  p->number = i;
	  p->next   = q;	  
	  p->back   = (node *) NULL;
	  p->hash   = 0;       
	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;  

  tr->ntips       = 0;
  tr->nextnode    = 0;
 
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;

  tr->bitVectors = (unsigned int **)NULL;

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;
  
  tr->nameHash = initStringHashTable(10 * tr->mxtips);

  tr->partitionData = (pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);

  return TRUE;
}


int partComparePtr(const void *p1, const void *p2)
{
  pInfo *p1p = *((pInfo**)p1),
    *p2p = *((pInfo**)p2); 

  return p1p->width - p2p->width ; 
}


/* TODO consolidate with partComparePtr */
static int partCompare(const void *p1, const void *p2)
{
 partitionType 
   *rc1 = (partitionType *)p1,
   *rc2 = (partitionType *)p2;

 int 
   i = rc1->partitionLength,
   j = rc2->partitionLength;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}



static void multiprocessorScheduling(tree *tr)
{
  int 
    s,
    model,
    modelStates[2] = {4, 20},
    numberOfPartitions[2] = {0 , 0},
    arrayLength = sizeof(modelStates) / sizeof(int);
  
    /* check that we have not addedd any new models for data types with a different number of states
       and forgot to update modelStates */
    
    tr->partitionAssignment = (int *)malloc(tr->NumberOfModels * sizeof(int));
    
  for(model = 0; model < tr->NumberOfModels; model++)
    {        
      boolean 
	exists = FALSE;

      for(s = 0; s < arrayLength; s++)
	{
	  exists = exists || (tr->partitionData[model].states == modelStates[s]);
	  if(tr->partitionData[model].states == modelStates[s])
	    numberOfPartitions[s] += 1;
	}

      assert(exists);
    }

  if(ABS_ID(tr) == 0 )
    printBothOpen(tr,"\nMulti-processor partition data distribution enabled (-Q option)\n");

  for(s = 0; s < arrayLength; s++)
    {
      if(numberOfPartitions[s] > 0)
	{
	  size_t   
	    checkSum = 0,
	    sum = 0;
	  
	  int    
	    i,
	    k,
	    n = ABS_NUM_RANK,
	    p = numberOfPartitions[s],    
	    *assignments = (int *)calloc(n, sizeof(int));  
	  
	  partitionType 
	    *pt = (partitionType *)malloc(sizeof(partitionType) * p);
	  
	  
	  for(i = 0, k = 0; i < tr->NumberOfModels; i++)
	    {
	      if(tr->partitionData[i].states == modelStates[s])
		{
		  pt[k].partitionNumber = i;
		  pt[k].partitionLength = tr->partitionData[i].upper - tr->partitionData[i].lower;
		  sum += (size_t)pt[k].partitionLength;
		  k++;
		}
	    }
	  
	  assert(k == p);
	  
	  qsort(pt, p, sizeof(partitionType), partCompare);    
	  
	  for(i = 0; i < p; i++)
	    {
	      int 
		k, 
		min = INT_MAX,
		minIndex = -1;
	      
	      for(k = 0; k < n; k++)	
		if(assignments[k] < min)
		  {
		    min = assignments[k];
		    minIndex = k;
		  }
	      
	      assert(minIndex >= 0);
	      
	      assignments[minIndex] +=  pt[i].partitionLength;
	      assert(pt[i].partitionNumber >= 0 && pt[i].partitionNumber < tr->NumberOfModels);
	      tr->partitionAssignment[pt[i].partitionNumber] = minIndex;
	    }

	  if(ABS_ID(tr) == 0)
	    {
	      for(i = 0; i < n; i++)	       
		printBothOpen(tr,"Process %d has %d sites for %d state model \n", i, assignments[i], modelStates[s]); 
	      
	      printBothOpen(tr,"\n");
	    }

	  for(i = 0; i < n; i++)
	    checkSum += (size_t)assignments[i];
	  
	  assert(sum == checkSum);
	  
	  free(assignments);
	  free(pt);
	}
    }
}

static void computeFraction(tree *tr)
{
  int
    model;

  size_t 
    i;

  int assigned = 0; 
  int unassigned = 0; 
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      size_t 
	width = 0;

      for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	{
	  if(i % ABS_NUM_RANK == ABS_ID(tr))
	    { 
	      width++;
	      assigned++; 
	    }
	  else 
	    {
	      unassigned++; 
	    }
	}

      tr->partitionData[model].width = width;
    }

  printf("assigned %d, unassigned %d\n", assigned, unassigned); 
}


static void computeFractionMany(tree *tr)
{
  int
    sites = 0;

  int   
    model;

  assert(tr->manyPartitions);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(isThisMyPartition(tr,model))
	{	 
	  tr->partitionData[model].width = tr->partitionData[model].upper - tr->partitionData[model].lower;
	  sites += tr->partitionData[model].width;
	}
      else       	  
	tr->partitionData[model].width = 0;       
    }  
}


static int iterated_bitcount(unsigned int n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}


static void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile)
{  
  size_t
    bytes_read;
  
  bytes_read = fread(ptr, size, nmemb, byteFile);

  assert(bytes_read == nmemb);
}

/*static char bits_in_16bits [0x1u << 16];*/

static void compute_bits_in_16bits(char *bits_in_16bits)
{
    unsigned int i;    
    
    /* size is 65536 */

    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);       

    return ;
}

unsigned int precomputed16_bitcount (unsigned int n, char *bits_in_16bits)
{
  /* works only for 32-bit int*/
    
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}



void initializePartitions(tree *tr, FILE *byteFile)
{ 
  size_t
    i,
    j,
    width,
    model,
    /* offset, */
    numRank = ABS_NUM_RANK, 
    myGlobalRank = ABS_ID(tr),
    countOffset,
    myLength = 0;

  int
    maxCategories;

  unsigned char 
    *y;

  compute_bits_in_16bits(tr->bits_in_16bits);

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    tr->partitionData[model].width        = 0;

  if(tr->manyPartitions)
    {
      multiprocessorScheduling(tr);  
      computeFractionMany(tr);
    }
  else
    computeFraction(tr);
  	   
  maxCategories = tr->maxCategories;

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {                       
      const partitionLengths 
	*pl = getPartitionLengths(&(tr->partitionData[model])); 

      width = tr->partitionData[model].width;

      tr->partitionData[model].wr = (double *)malloc(sizeof(double) * width);
      tr->partitionData[model].wr2 = (double *)malloc(sizeof(double) * width);     

     	
      /* 
	 globalScaler needs to be 2 * tr->mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> calloc
       */

      tr->partitionData[model].globalScaler    = (unsigned int *)calloc(2 * tr->mxtips, sizeof(unsigned int));  	         

      tr->partitionData[model].left              = (double *)malloc_aligned(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].right             = (double *)malloc_aligned(pl->rightLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].EIGN              = (double*)malloc(pl->eignLength * sizeof(double));
      tr->partitionData[model].EV                = (double*)malloc_aligned(pl->evLength * sizeof(double));
      tr->partitionData[model].EI                = (double*)malloc(pl->eiLength * sizeof(double));
      
      tr->partitionData[model].substRates        = (double *)malloc(pl->substRatesLength * sizeof(double));
      tr->partitionData[model].frequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].empiricalFrequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].tipVector         = (double *)malloc_aligned(pl->tipVectorLength * sizeof(double));
      tr->partitionData[model].symmetryVector    = (int *)malloc(pl->symmetryVectorLength  * sizeof(int));
      tr->partitionData[model].frequencyGrouping = (int *)malloc(pl->frequencyGroupingLength  * sizeof(int));
      
      tr->partitionData[model].perSiteRates      = (double *)malloc(sizeof(double) * tr->maxCategories);
            
      tr->partitionData[model].nonGTR = FALSE;            

      tr->partitionData[model].gammaRates = (double*)malloc(sizeof(double) * 4);
      tr->partitionData[model].yVector = (unsigned char **)malloc(sizeof(unsigned char*) * (tr->mxtips + 1));

      
      tr->partitionData[model].xVector = (double **)malloc(sizeof(double*) * tr->mxtips);   
      	
      for(j = 0; j < (size_t)tr->mxtips; j++)	        	  	  	  	 
	  tr->partitionData[model].xVector[j]   = (double*)NULL;   

      tr->partitionData[model].xSpaceVector = (size_t *)calloc(tr->mxtips, sizeof(size_t));  

      tr->partitionData[model].sumBuffer = (double *)malloc_aligned(width *
									   (size_t)(tr->partitionData[model].states) *
									   discreteRateCategories(tr->rateHetModel) *
									   sizeof(double));
	    
      tr->partitionData[model].wgt = (int *)malloc_aligned(width * sizeof(int));	  

      /* rateCategory must be assigned using calloc() at start up there is only one rate category 0 for all sites */
      
      /* :TODO: reduce memory consumption? => has a significant influence on vectorized stuff later, careful!   */
      tr->partitionData[model].rateCategory = (int *)calloc(width, sizeof(int));
      /* tr->partitionData[model].rateCategory = (int *)calloc(width, sizeof(unsigned char)); */

      if(width > 0 && tr->saveMemory)
	{
	  tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
	    
	  tr->partitionData[model].gapVector = (unsigned int*)calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));	  	    	  	  
	    
	  tr->partitionData[model].gapColumn = (double *)malloc_aligned(((size_t)tr->mxtips) *								      
									       ((size_t)(tr->partitionData[model].states)) *
									       discreteRateCategories(tr->rateHetModel) * sizeof(double));
	}
      else
	{
	   tr->partitionData[model].gapVectorLength = 0;
	    
	   tr->partitionData[model].gapVector = (unsigned int*)NULL; 	  	    	   
	    
	   tr->partitionData[model].gapColumn = (double*)NULL;	    	    	   
	}              
    }

        
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    myLength += tr->partitionData[model].width;         
   
  /* assign local memory for storing sequence data */

  tr->y_ptr = (unsigned char *)malloc(myLength * (size_t)(tr->mxtips) * sizeof(unsigned char));
  assert(tr->y_ptr != NULL);
   
  for(i = 0; i < (size_t)tr->mxtips; i++)
    {
      for(model = 0,countOffset = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  tr->partitionData[model].yVector[i+1]   = &tr->y_ptr[i * myLength + countOffset];
	  countOffset +=  tr->partitionData[model].width;
	}
      assert(countOffset == myLength);
    }

  /* figure in data */

  if(tr->manyPartitions)
    {
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  if(isThisMyPartition(tr, model))
	    {
	      width = tr->partitionData[model].upper - tr->partitionData[model].lower;	     
	      
	      memcpy(&(tr->partitionData[model].wgt[0]), &(tr->aliaswgt[tr->partitionData[model].lower]), sizeof(int) * width);
	    }
	}
    }
  else
    {
      size_t 	   
	globalCounter, 
	r, 
	localCounter;
      
      for(model = 0, globalCounter = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  for(localCounter = 0, r = (size_t)tr->partitionData[model].lower;  r < (size_t)tr->partitionData[model].upper; r++)
	    {
	      if(r % numRank == myGlobalRank)
		{
		  tr->partitionData[model].wgt[localCounter] = tr->aliaswgt[globalCounter];	      	     		 		  					     
		  
		  localCounter++;
		}
	      globalCounter++;
	    }
	  assert(localCounter == tr->partitionData[model].width);
	}   
      assert(globalCounter == tr->originalCrunchedLength);
    }
   
  y = (unsigned char *)malloc(sizeof(unsigned char) * tr->originalCrunchedLength);

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength), byteFile);
	
      if(tr->manyPartitions)
	{
	  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	    {
	      if(isThisMyPartition(tr, model))
		{
		  memcpy(tr->partitionData[model].yVector[i], &(y[tr->partitionData[model].lower]), sizeof(unsigned char) * tr->partitionData[model].width);					    
		  assert(tr->partitionData[model].width == tr->partitionData[model].upper - tr->partitionData[model].lower);
		}
	      else
		assert(tr->partitionData[model].width == 0);
	    }
	}
      else
	{
	  size_t 	  
	    globalCounter, 
	    r, 
	    localCounter;

	  for(model = 0, globalCounter = 0; model < (size_t)tr->NumberOfModels; model++)
	    {
	      for(localCounter = 0, r = (size_t)tr->partitionData[model].lower;  r < (size_t)tr->partitionData[model].upper; r++)
		{
		  if(r % numRank ==  myGlobalRank )
		    {		      
		      tr->partitionData[model].yVector[i][localCounter] = y[globalCounter]; 	     
		      
		      localCounter++;
		    }
		  globalCounter++;
		}
	      
	      assert(localCounter == tr->partitionData[model].width);
	    }

	  assert(globalCounter == tr->originalCrunchedLength);
	}
    }

  free(y);
    
  /* initialize gap bit vectors at tips when memory saving option is enabled */
  
  if(tr->saveMemory)
    {
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  int        
	    undetermined = getUndetermined(tr->partitionData[model].dataType);
	  	 
	  width =  tr->partitionData[model].width;
	    
	  if(width > 0)
	    {	   	    	      	    	     
	      for(j = 1; j <= (size_t)(tr->mxtips); j++)
		for(i = 0; i < width; i++)
		  if(tr->partitionData[model].yVector[j][i] == undetermined)
		    tr->partitionData[model].gapVector[tr->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];	    
	    }     
	}
    }
}




const partitionLengths*  getPartitionLengths(pInfo *p)
{
  /* :TODO: this can easily be re-hacked  */
  /* assert(0);  */

  int
    dataType  = p->dataType
    /* ,states    = p->states */
    /* ,tipLength = p->maxTipStates */
;

  /* assert(states != -1 && tipLength != -1); */

  /* assert(MIN_MODEL < dataType && dataType < MAX_MODEL); */

  /* pLength.leftLength = pLength.rightLength = states * states; */
  /* pLength.eignLength = states; */
  /* pLength.evLength   = states * states; */
  /* pLength.eiLength   = states * states; */
  /* pLength.substRatesLength = (states * states - states) / 2; */
  /* pLength.frequenciesLength = states; */
  /* pLength.tipVectorLength   = tipLength * states; */
  /* pLength.symmetryVectorLength = (states * states - states) / 2; */
  /* pLength.frequencyGroupingLength = states; */
  /* pLength.nonGTR = FALSE; */

  return (&pLengths[dataType]); 
}



void initializeTree(tree *tr, analdef *adef)
{ 
  size_t 
    i,
    model;
  
  FILE 
    *byteFile = fopen(byteFileName, "rb");

  double 
    **empiricalFrequencies;	 

  myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
  myBinFread(&(tr->originalCrunchedLength), sizeof(size_t), 1, byteFile);
  myBinFread(&(tr->NumberOfModels),         sizeof(int), 1, byteFile);
  myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);
   
  empiricalFrequencies = (double **)malloc(sizeof(double *) * tr->NumberOfModels);
  
  if(adef->perGeneBranchLengths)
    tr->numBranches = tr->NumberOfModels;
  else
    tr->numBranches = 1;
  
  /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
     let's not worry about this right now, because it is indeed ExaML-specific */
  
 
    
  tr->aliaswgt                   = (int *)malloc(tr->originalCrunchedLength * sizeof(int));
  myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);	       
  
  tr->rateCategory    = (int *)    calloc(tr->originalCrunchedLength, sizeof(int));	  
  tr->wr              = (double *) malloc(tr->originalCrunchedLength * sizeof(double)); 
  tr->wr2             = (double *) malloc(tr->originalCrunchedLength * sizeof(double)); 
  tr->patrat          = (double*)  malloc(tr->originalCrunchedLength * sizeof(double));
  tr->patratStored    = (double*)  malloc(tr->originalCrunchedLength * sizeof(double)); 
  tr->lhs             = (double*)  malloc(tr->originalCrunchedLength * sizeof(double)); 
  
  tr->executeModel   = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);
  
  for(i = 0; i < (size_t)tr->NumberOfModels; i++)
    tr->executeModel[i] = TRUE;
   
  setupTree(tr); 

  if(tr->searchConvergenceCriterion && ABS_ID(tr) == 0)
    { 
      tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
      tr->h = initHashTable(tr->mxtips * 4);     
    }
  
  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      int 
	len;
      
      myBinFread(&len, sizeof(int), 1, byteFile);
      tr->nameList[i] = (char*)malloc(sizeof(char) * len);
      myBinFread(tr->nameList[i], sizeof(char), len, byteFile);
      addword(tr->nameList[i], tr->nameHash, i);        
    }  
 
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {      
      int 
	len;
      
      pInfo 
	*p = &(tr->partitionData[model]);	   
      
      myBinFread(&(p->states),             sizeof(int), 1, byteFile);
      myBinFread(&(p->maxTipStates),       sizeof(int), 1, byteFile);
      myBinFread(&(p->lower),              sizeof(size_t), 1, byteFile);
      myBinFread(&(p->upper),              sizeof(size_t), 1, byteFile);
      myBinFread(&(p->width),              sizeof(size_t), 1, byteFile);
      myBinFread(&(p->dataType),           sizeof(int), 1, byteFile);
      myBinFread(&(p->protModels),         sizeof(int), 1, byteFile);
      myBinFread(&(p->autoProtModels),     sizeof(int), 1, byteFile);
      myBinFread(&(p->protFreqs),          sizeof(int), 1, byteFile);
      myBinFread(&(p->nonGTR),             sizeof(boolean), 1, byteFile);
      myBinFread(&(p->numberOfCategories), sizeof(int), 1, byteFile); 
      /* later on if adding secondary structure data
	 
	 int    *symmetryVector;
	 int    *frequencyGrouping;
      */
      
      myBinFread(&len, sizeof(int), 1, byteFile);
      p->partitionName = (char*)malloc(sizeof(char) * len);
      myBinFread(p->partitionName, sizeof(char), len, byteFile);
      
      empiricalFrequencies[model] = (double *)malloc(sizeof(double) * p->states);
      myBinFread(empiricalFrequencies[model], sizeof(double), p->states, byteFile);	   
    }
  
  initializePartitions(tr, byteFile);

  fclose(byteFile);


  initModel(tr, empiricalFrequencies); 

  tr->reductionBuffer = calloc(2 * tr->numBranches, sizeof(double)); 
 
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    free(empiricalFrequencies[model]);

  free(empiricalFrequencies);
}
