/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#include "axml.h"
#include "tree.h" 

#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
/*
  special bug fix, enforces denormalized numbers to be flushed to zero,
  without this program is a tiny bit faster though.
  #include <emmintrin.h> 
  #define MM_DAZ_MASK    0x0040
  #define MM_DAZ_ON    0x0040
  #define MM_DAZ_OFF    0x0000
*/
#endif



#define INCLUDE_DEFINITION
#include "globalVariables.h"
#undef INCLUDE_DEFINITION


/***************** UTILITY FUNCTIONS **************************/

void storeExecuteMaskInTraversalDescriptor(tree *tr)
{
   int model;
      
   for(model = 0; model < tr->NumberOfModels; model++)
     tr->td[0].executeModel[model] = tr->executeModel[model];
}

void storeValuesInTraversalDescriptor(tree *tr, double *value)
{
   int model;
      
   for(model = 0; model < tr->NumberOfModels; model++)
     tr->td[0].parameterValues[model] = value[model];
}


void *malloc_aligned(size_t size) 
{
  void 
    *ptr = (void *)NULL;
 
  int 
    res;
  

#if defined (__APPLE__)
  /* 
     presumably malloc on MACs always returns 
     a 16-byte aligned pointer
  */

  ptr = malloc(size);
  
  if(ptr == (void*)NULL) 
   assert(0);
  
#ifdef __AVX
  assert(0);
#endif


#else
  res = posix_memalign( &ptr, BYTE_ALIGNMENT, size );

  if(res != 0) 
    assert(0);
#endif 
   
  return ptr;
}






static void printBoth(FILE *f, const char* format, ... )
{
  if(mpiState.rank == 0)
    {
      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
    }
}

void printBothOpen(const char* format, ... )
{
  if(mpiState.rank == 0)
    {
      FILE *f = myfopen(infoFileName, "ab");
      
      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
      
      fclose(f);
    }
}




boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}


int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}

int getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}



char getInverseMeaning(int dataType, unsigned char state)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return  pLengths[dataType].inverseMeaning[state];
}

partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
    states    = p->states,
    tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  pLength.leftLength = pLength.rightLength = states * states;
  pLength.eignLength = states;
  pLength.evLength   = states * states;
  pLength.eiLength   = states * states;
  pLength.substRatesLength = (states * states - states) / 2;
  pLength.frequenciesLength = states;
  pLength.tipVectorLength   = tipLength * states;
  pLength.symmetryVectorLength = (states * states - states) / 2;
  pLength.frequencyGroupingLength = states;
  pLength.nonGTR = FALSE;

  return (&pLengths[dataType]); 
}










size_t discreteRateCategories(int rateHetModel)
{
  size_t 
    result;

  switch(rateHetModel)
    {
    case CAT:
      result = 1;
      break;
    case GAMMA:
      result = 4;
      break;
    default:
      assert(0);
    }

  return result;
}



double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

static int filexists(char *filename)
{
  FILE 
    *fp;
  
  int 
    res;
  
  fp = fopen(filename,"rb");

  if(fp)
    {
      res = 1;
      fclose(fp);
    }
  else
    res = 0;

  return res;
}


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {
      if(fp)
	return fp;
      else
	{
	  if(mpiState.rank == 0)
	    printf("The file %s you want to open for reading does not exist, exiting ...\n", path);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	  if(mpiState.rank == 0)
	    printf("The file %s ExaML wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
		   path, mode);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }


}





/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/


boolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}









void getxnode (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }

  assert(p->x);
}





void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;
}


/***********************reading and initializing input ******************/







boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}



static void initAdef(analdef *adef)
{   
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = 10;
  adef->bestTrav               = 10;
  adef->initialSet             = FALSE; 
  adef->mode                   = BIG_RAPID_MODE; 
  adef->likelihoodEpsilon      = 0.1;
 
  adef->permuteTreeoptimize    = FALSE; 
  adef->perGeneBranchLengths   = FALSE;  
 
  adef->useCheckpoint          = FALSE;
   
#ifdef _BAYESIAN 
  adef->bayesian               = FALSE;
#endif

}




static int modelExists(char *model, tree *tr)
{  
   if(strcmp(model, "PSR\0") == 0)
    {
      tr->rateHetModel = CAT;
      return 1;
    }

  if(strcmp(model, "GAMMA\0") == 0)
    {
      tr->rateHetModel = GAMMA;
      return 1;
    }

  
  return 0;
}



static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
	return -1;
    }
  else
    {
      if(strcmp(argv[*optind], "--") == 0)
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0)
    {
      printf(": illegal option -- %c \n", c);
      if(argv[*optind][++sp] == '\0')
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      return('?');
    }
  if(*++cp == ':')
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc)
	    {
	      printf(": option requires an argument -- %c\n", c);
	      sp = 1;
	      return('?');
	    }
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    }
  else
    {
      if(argv[*optind][++sp] == '\0')
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
  return(c);
  }




/*********************************** *********************************************************/


static void printVersionInfo(void)
{
  if(mpiState.rank == 0)
    printf("\n\nThis is %s version %s released by Alexandros Stamatakis on %s.\n\n",  programName, programVersion, programDate); 
}

static void printMinusFUsage(void)
{
  printf("\n");
 

  printf("              \"-f d\": new rapid hill-climbing \n");
  printf("                      DEFAULT: ON\n");

  printf("\n");

  printf("              \"-f o\": old and slower rapid hill-climbing without heuristic cutoff\n");
  
  printf("\n");

  printf("              DEFAULT for \"-f\": new rapid hill climbing\n");

  printf("\n");
}


static void printREADME(void)
{
  if(mpiState.rank == 0)
    {
      printVersionInfo();
      printf("\n");  
      printf("\nTo report bugs use the RAxML google group\n");
      printf("Please send me all input files, the exact invocation, details of the HW and operating system,\n");
      printf("as well as all error messages printed to screen.\n\n\n");
      
      printf("examl|examl-AVX\n");
      printf("      -s binarySequenceFileName\n");
      printf("      -n outputFileNames\n");
      printf("      -m rateHeterogeneityModel\n");
      printf("      -t userStartingTree|-R binaryCheckpointFile\n");
      printf("      [-a]\n");
      printf("      [-B numberOfMLtreesToSave]\n"); 
      printf("      [-c numberOfCategories]\n");
      printf("      [-D]\n");
      printf("      [-e likelihoodEpsilon] \n");
      printf("      [-f d|o]\n");    
      printf("      [-h] \n");
      printf("      [-i initialRearrangementSetting] \n");
      printf("      [-M]\n");
      printf("      [-Q]\n");
      printf("      [-S]\n");
      printf("      [-v]\n"); 
      printf("      [-w outputDirectory] \n"); 
      printf("\n");  
      printf("      -a      use the median for the discrete approximation of the GAMMA model of rate heterogeneity\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");
      printf("\n");
      printf("      -B      specify the number of best ML trees to save and print to file\n");
      printf("\n");
      printf("      -c      Specify number of distinct rate catgories for ExaML when modelOfEvolution\n");
      printf("              is set to GTRPSR\n");
      printf("              Individual per-site rates are categorized into numberOfCategories rate \n");
      printf("              categories to accelerate computations. \n");
      printf("\n");
      printf("              DEFAULT: 25\n");
      printf("\n");
      printf("      -D      ML search convergence criterion. This will break off ML searches if the relative \n");
      printf("              Robinson-Foulds distance between the trees obtained from two consecutive lazy SPR cycles\n");
      printf("              is smaller or equal to 1%s. Usage recommended for very large datasets in terms of taxa.\n", "%");
      printf("              On trees with more than 500 taxa this will yield execution time improvements of approximately 50%s\n",  "%");
      printf("              While yielding only slightly worse trees.\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");    
      printf("\n");
      printf("      -e      set model optimization precision in log likelihood units for final\n");
      printf("              optimization of model parameters\n");
      printf("\n");
      printf("              DEFAULT: 0.1 \n"); 
      printf("\n");
      printf("      -f      select algorithm:\n");
      
      printMinusFUsage();
 
      printf("\n");
      printf("      -h      Display this help message.\n");
      printf("\n");  
      printf("      -i      Initial rearrangement setting for the subsequent application of topological \n");
      printf("              changes phase\n");
      printf("\n");
      printf("      -m      Model of rate heterogeneity\n");
      printf("\n"); 
      printf("              select \"-m PSR\" for the per-site rate category model (this used to be called CAT in RAxML)\n");
      printf("              select \"-m GAMMA\" for the gamma model of rate heterogeneity with 4 discrete rates\n");
      printf("\n");
      printf("      -M      Switch on estimation of individual per-partition branch lengths. Only has effect when used in combination with \"-q\"\n");
      printf("              Branch lengths for individual partitions will be printed to separate files\n");
      printf("              A weighted average of the branch lengths is computed by using the respective partition lengths\n");
      printf("\n");
      printf("              DEFAULT: OFF\n");
      printf("\n");
      printf("      -n      Specifies the name of the output file.\n"); 
      printf("\n");
      printf("      -Q      Enable alternative data/load distribution algorithm for datasets with many partitions\n");
      printf("              In particular under PSR this can lead to parallel performance improvements of up to factor 10!\n");
      printf("\n");
      printf("      -R      read in a binary checkpoint file called ExaML_binaryCheckpoint.RUN_ID_number\n");
      printf("\n");
      printf("      -s      Specify the name of the BINARY alignment data file generated by the parser component\n");
      printf("\n");
      printf("      -S      turn on memory saving option for gappy multi-gene alignments. For large and gappy datasets specify -S to save memory\n");
      printf("              This will produce slightly different likelihood values, may be a bit slower but can reduce memory consumption\n");
      printf("              from 70GB to 19GB on very large and gappy datasets\n");
      printf("\n");
      printf("      -t      Specify a user starting tree file name in Newick format\n");
      printf("\n");
      printf("      -v      Display version information\n");
      printf("\n");
      printf("      -w      FULL (!) path to the directory into which ExaML shall write its output files\n");
      printf("\n");
      printf("              DEFAULT: current directory\n");  
      printf("\n\n\n\n");
    }
}




static void analyzeRunId(char id[128])
{
  int i = 0;

  while(id[i] != '\0')
    {    
      if(i >= 128)
	{
	  printf("Error: run id after \"-n\" is too long, it has %d characters please use a shorter one\n", i);
	  assert(0);
	}
      
      if(id[i] == '/')
	{
	  printf("Error character %c not allowed in run ID\n", id[i]);
	  assert(0);
	}


      i++;
    }

  if(i == 0)
    {
      printf("Error: please provide a string for the run id after \"-n\" \n");
      assert(0);
    }

}

static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  boolean
    bad_opt    =FALSE,
    resultDirSet = FALSE;

  char
    resultDir[1024] = "",          
    *optarg,
    model[1024] = "",       
    modelChar;

  double 
    likelihoodEpsilon;
  
  int     
    optind = 1,        
    c,
    nameSet = 0,
    treeSet = 0,   
    modelSet = 0, 
    byteFileSet = 0;

#ifdef _USE_PTHREADS
  NumberOfThreads = 0;
#endif


  /*********** tr inits **************/ 
 
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA;
 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;

  tr->manyPartitions = FALSE;

  tr->categories             = 25;

  tr->grouped = FALSE;
  tr->constrained = FALSE;

  tr->gapyness               = 0.0; 
  tr->saveBestTrees          = 0;

  tr->useMedian = FALSE;

  /********* tr inits end*************/


  

  while(!bad_opt && ((c = mygetopt(argc,argv,"R:B:e:c:f:i:m:t:T:w:n:s:vhMSDQa", &optind, &optarg))!=-1))
    {
    switch(c)
      {    
      case 'a':
	tr->useMedian = TRUE;
	break;
      case 'B':
	sscanf(optarg,"%d", &(tr->saveBestTrees));
	if(tr->saveBestTrees < 0)
	  {
	    printf("Number of best trees to save must be greater than 0!\n");
	    errorExit(-1);	 
	  }
	break;       
      case 'Q':
	tr->manyPartitions = TRUE; 
	break;
      case 's':		 	
	strcpy(byteFileName, optarg);	 	
	byteFileSet = TRUE;
	/*printf("%s \n", byteFileName);*/
	break;      
      case 'S':
	tr->saveMemory = TRUE;
	break;
      case 'D':
	tr->searchConvergenceCriterion = TRUE;	
	break;
      case 'R':
	adef->useCheckpoint = TRUE;
	strcpy(binaryCheckpointInputName, optarg);
	break;          
      case 'M':
	adef->perGeneBranchLengths = TRUE;
	break;                                 
      case 'e':
	sscanf(optarg,"%lf", &likelihoodEpsilon);
	adef->likelihoodEpsilon = likelihoodEpsilon;
	break;    
      
      case 'v':
	printVersionInfo();
	errorExit(0);
      
      case 'h':
	printREADME();
	errorExit(0);     
      case 'c':
	sscanf(optarg, "%d", &tr->categories);
	break;     
      case 'f':
	sscanf(optarg, "%c", &modelChar);
	switch(modelChar)
	  {	 
	  case 'd':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;	  
	  case 'o':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = FALSE;
	    break;	    	  	  	     
	  default:
	    {
	      if(mpiState.rank == 0)
		{
		  printf("Error select one of the following algorithms via -f :\n");
		  printMinusFUsage();
		}
	      errorExit(-1);
	    }
	  }
	break;
      case 'i':
	sscanf(optarg, "%d", &adef->initial);
	adef->initialSet = TRUE;
	break;
      case 'n':
        strcpy(run_id,optarg);
	analyzeRunId(run_id);
	nameSet = 1;
        break;
      case 'w':
        strcpy(resultDir, optarg);
	resultDirSet = TRUE;
        break;
      case 'T':
#ifdef _USE_PTHREADS
	sscanf(optarg,"%d", &NumberOfThreads);
#else
	if(mpiState.rank == 0)
	  {
	    printf("Option -T does not have any effect with the sequential or parallel MPI version.\n");
	    printf("It is used to specify the number of threads for the Pthreads-based parallelization\n");
	  }	
#endif
	break; 
      case 't':
	strcpy(tree_file, optarg);       
	treeSet = 1;       
	break;     
      case 'm':
	strcpy(model,optarg);
	if(modelExists(model, tr) == 0)
	  {
	    if(mpiState.rank == 0)
	      {
		printf("Rate heterogeneity Model %s does not exist\n\n", model);               
		printf("For per site rates (called CAT in previous versions) use: PSR\n");	
		printf("For GAMMA use: GAMMA\n");		
	      }
	    errorExit(-1);
	  }
	else
	  modelSet = 1;
	break;     
      default:
	errorExit(-1);
      }
    }


  if(!byteFileSet)
    {
      if(mpiState.rank == 0)
	printf("\nError, you must specify a binary format data file with the \"-s\" option\n");
      errorExit(-1);
    }

  if(!modelSet)
    {
      if(mpiState.rank == 0)
	printf("\nError, you must specify a model of rate heterogeneity with the \"-m\" option\n");
      errorExit(-1);
    }

  if(!nameSet)
    {
      if(mpiState.rank == 0)
	printf("\nError: please specify a name for this run with -n\n");
      errorExit(-1);
    }

  if(!treeSet && !adef->useCheckpoint)
    {
      if(mpiState.rank == 0)
	{
	  printf("\nError: please either specify a starting tree for this run with -t\n");
	  printf("or re-start the run from a checkpoint with -R\n");
	}
      
      errorExit(-1);
    }
  
   {

    const 
      char *separator = "/";

    if(resultDirSet)
      {
	char 
	  dir[1024] = "";
	

	if(resultDir[0] != separator[0])
	  strcat(dir, separator);
	
	strcat(dir, resultDir);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	strcpy(workdir, dir);
      }
    else
      {
	char 
	  dir[1024] = "",
	  *result = getcwd(dir, sizeof(dir));
	
	assert(result != (char*)NULL);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	
	strcpy(workdir, dir);		
      }
   }

  return;
}




void errorExit(int e)
{
  MPI_Finalize();

  exit(e);
}



static void makeFileNames(void)
{
  int 
    infoFileExists = 0;
    
  strcpy(resultFileName,       workdir);
  strcpy(logFileName,          workdir);  
  strcpy(infoFileName,         workdir);
  strcpy(binaryCheckpointName, workdir);
   
  strcat(resultFileName,       "ExaML_result.");
  strcat(logFileName,          "ExaML_log.");  
  strcat(infoFileName,         "ExaML_info.");
  strcat(binaryCheckpointName, "ExaML_binaryCheckpoint.");
  
  strcat(resultFileName,       run_id);
  strcat(logFileName,          run_id);  
  strcat(infoFileName,         run_id); 
  strcat(binaryCheckpointName, run_id);

  infoFileExists = filexists(infoFileName);

  if(infoFileExists)
    {
      if(mpiState.rank == 0)
	{
	  printf("ExaML output files with the run ID <%s> already exist \n", run_id);
	  printf("in directory %s ...... exiting\n", workdir);
	}

#ifndef _NOT_PRODUCTIVE
      errorExit(-1);	
#endif
    }
}




 




/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/


static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{
  if(mpiState.rank == 0)
    {
      int i, model;
      FILE *infoFile = myfopen(infoFileName, "ab");
      char modelType[128];

      
      if(tr->useMedian)
	strcpy(modelType, "GAMMA with Median");
      else
	strcpy(modelType, "GAMMA");   
     
      printBoth(infoFile, "\n\nThis is %s version %s released by Alexandros Stamatakis in %s.\n\n",  programName, programVersion, programDate);
     
      
      
     
      printBoth(infoFile, "\nAlignment has %d distinct alignment patterns\n\n",  tr->originalCrunchedLength);
      
     
      
      printBoth(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * tr->gapyness, "%");
      

      switch(adef->mode)
	{	
	case  BIG_RAPID_MODE:	 
	  printBoth(infoFile, "\nExaML rapid hill-climbing mode\n\n");
	  break;	
	default:
	  assert(0);
	}

     
	  
      if(adef->perGeneBranchLengths)
	printBoth(infoFile, "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", tr->NumberOfModels);
      else
	printBoth(infoFile, "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", tr->NumberOfModels);	
	

      
     


      
      
      printBoth(infoFile, "All free model parameters will be estimated by ExaML\n");
      
     
	
      if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
	printBoth(infoFile, "%s model of rate heteorgeneity, ML estimate of alpha-parameter\n\n", modelType);
      else
	{
	  printBoth(infoFile, "ML estimate of %d per site rate categories\n\n", tr->categories);
	  /*
	    if(adef->mode != CLASSIFY_ML)
	    printBoth(infoFile, "Likelihood of final tree will be evaluated and optimized under %s\n\n", modelType);
	  */
	}
      
      /*
	if(adef->mode != CLASSIFY_ML)
	printBoth(infoFile, "%s Model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n",
	modelType, adef->likelihoodEpsilon);
      */
    
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  printBoth(infoFile, "Partition: %d\n", model);
	  printBoth(infoFile, "Alignment Patterns: %d\n", tr->partitionData[model].upper - tr->partitionData[model].lower);
	  printBoth(infoFile, "Name: %s\n", tr->partitionData[model].partitionName);
	  
	  switch(tr->partitionData[model].dataType)
	    {
	    case DNA_DATA:
	      printBoth(infoFile, "DataType: DNA\n");	     
	      printBoth(infoFile, "Substitution Matrix: GTR\n");
	      break;
	    case AA_DATA:
	      assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);
	      printBoth(infoFile, "DataType: AA\n");	      
	      printBoth(infoFile, "Substitution Matrix: %s\n", protModels[tr->partitionData[model].protModels]);
	      printBoth(infoFile, "%s Base Frequencies:\n", (tr->partitionData[model].protFreqs == 1)?"Empirical":"Fixed");	     
	      break;
	    case BINARY_DATA:
	      printBoth(infoFile, "DataType: BINARY/MORPHOLOGICAL\n");	      
	      printBoth(infoFile, "Substitution Matrix: Uncorrected\n");
	      break;
	    case SECONDARY_DATA:
	      printBoth(infoFile, "DataType: SECONDARY STRUCTURE\n");	     
	      printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	      break;
	    case SECONDARY_DATA_6:
	      printBoth(infoFile, "DataType: SECONDARY STRUCTURE 6 STATE\n");	     
	      printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	      break;
	    case SECONDARY_DATA_7:
	      printBoth(infoFile, "DataType: SECONDARY STRUCTURE 7 STATE\n");	      
	      printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	      break;
	    case GENERIC_32:
	      printBoth(infoFile, "DataType: Multi-State with %d distinct states in use (maximum 32)\n",tr->partitionData[model].states);		  
	      switch(tr->multiStateModel)
		{
		case ORDERED_MULTI_STATE:
		  printBoth(infoFile, "Substitution Matrix: Ordered Likelihood\n");
		  break;
		case MK_MULTI_STATE:
		  printBoth(infoFile, "Substitution Matrix: MK model\n");
		  break;
		case GTR_MULTI_STATE:
		  printBoth(infoFile, "Substitution Matrix: GTR\n");
		  break;
		default:
		  assert(0);
		}
	      break;
	    case GENERIC_64:
	      printBoth(infoFile, "DataType: Codon\n");		  
	      break;		
	    default:
	      assert(0);
	    }
	  printBoth(infoFile, "\n\n\n");
	}
      
      printBoth(infoFile, "\n");

      printBoth(infoFile, "ExaML was called as follows:\n\n");
      for(i = 0; i < argc; i++)
	printBoth(infoFile,"%s ", argv[i]);
      printBoth(infoFile,"\n\n\n");

      fclose(infoFile);
    }
}

void printResult(tree *tr, analdef *adef, boolean finalPrint)
{
  if(mpiState.rank == 0)
    {
      FILE *logFile;
      char temporaryFileName[1024] = "";
      
      strcpy(temporaryFileName, resultFileName);
      
      switch(adef->mode)
	{    
	case TREE_EVALUATION:
	  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, SUMMARIZE_LH, FALSE, FALSE);
	  
	  logFile = myfopen(temporaryFileName, "wb");
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);
	  
	  if(adef->perGeneBranchLengths)
	    printTreePerGene(tr, adef, temporaryFileName, "wb");
	  break;
	case BIG_RAPID_MODE:     
	  if(finalPrint)
	    {
	      switch(tr->rateHetModel)
		{
		case GAMMA:
		case GAMMA_I:
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint,
			      SUMMARIZE_LH, FALSE, FALSE);
		  
		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
		  
		  if(adef->perGeneBranchLengths)
		    printTreePerGene(tr, adef, temporaryFileName, "wb");
		  break;
		case CAT:
		  /*Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef,
		    NO_BRANCHES, FALSE, FALSE);*/
		  
		  
		  
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE,
			      TRUE, SUMMARIZE_LH, FALSE, FALSE);
		  
		  
		  
		  
		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
		  
		  break;
		default:
		  assert(0);
		}
	    }
	  else
	    {
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint,
			  NO_BRANCHES, FALSE, FALSE);
	      logFile = myfopen(temporaryFileName, "wb");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);
	    }    
	  break;
	default:
	  printf("FATAL ERROR call to printResult from undefined STATE %d\n", adef->mode);
	  exit(-1);
	  break;
	}
    }
}








void printLog(tree *tr)
{
  if(mpiState.rank == 0)
    {
      FILE *logFile;
      double t;
      
      t = gettime() - masterTime;
      
      logFile = myfopen(logFileName, "ab");
      
      /* printf("%f %1.40f\n", t, tr->likelihood); */

      fprintf(logFile, "%f %f\n", t, tr->likelihood);
      
      fclose(logFile);
    }
	     
}









void getDataTypeString(tree *tr, int model, char typeOfData[1024])
{
  switch(tr->partitionData[model].dataType)
    {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;
    case SECONDARY_DATA:
      strcpy(typeOfData,"SECONDARY 16 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_6:
      strcpy(typeOfData,"SECONDARY 6 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_7:
      strcpy(typeOfData,"SECONDARY 7 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
    default:
      assert(0);
    }
}




static void finalizeInfoFile(tree *tr, analdef *adef)
{
  if(mpiState.rank == 0)
    {
      double t;

      t = gettime() - masterTime;
      accumulatedTime = accumulatedTime + t;

      switch(adef->mode)
	{	
	case  BIG_RAPID_MODE:	 
	  printBothOpen("\n\nOverall Time for 1 Inference %f\n", t);
	  printBothOpen("\nOverall accumulated Time (in case of restarts): %f\n\n", accumulatedTime);
	  printBothOpen("Likelihood   : %f\n", tr->likelihood);
	  printBothOpen("\n\n");	  	  
	  printBothOpen("Final tree written to:                 %s\n", resultFileName);
	  printBothOpen("Execution Log File written to:         %s\n", logFileName);
	  printBothOpen("Execution information file written to: %s\n",infoFileName);	
	  break;
	default:
	  assert(0);
	}

	 
    }

}


/************************************************************************************/


/* 
	 tr->manyPartitions is set to TRUE if the user has indicated via -Q that there are substantially more partitions 
	 than threads/cores available. In that case we do not distribute sites from each partition in a cyclic fashion to the cores 
	 , but distribute entire partitions to cores. 
	 Achieving a good balance of alignment sites to cores boils down to the mult-processor scheduling problem known from theoretical comp. sci.
	 which si NP-complete.
	 We have implemented very simple "standard" heuristics for solving the multiprocessor scheduling problem that turn out to work very well
	 and are cheap to compute 
*/


static void examl_initMPI(int argc, char **argv)
{
  /* :TODO: catch errors in this function  */

  MPI_Init(&argc, &argv);
  mpiState.comm = MPI_COMM_WORLD; 
  MPI_Comm_rank(mpiState.comm, &mpiState.rank); 
  MPI_Comm_size(mpiState.comm, &mpiState.commSize);
  mpiState.mpiError = MPI_SUCCESS;
  mpiState.generation[PHASE_BRANCH_OPT] = 0; 
  mpiState.generation[PHASE_LNL_EVAL] = 0; 
  mpiState.generation[PHASE_RATE_OPT] = 0; 

#ifdef _USE_RTS 
  MPI_Comm_set_errhandler(mpiState.comm, MPI_ERRORS_RETURN);
  puts("Running with RTS activated"); 
#endif
  
  /* :TODO:necessary?  */
  MPI_Barrier(mpiState.comm);
}



int main (int argc, char *argv[])
{ 
  examl_initMPI(argc, argv);
  
  {
    tree  *tr = (tree*)malloc(sizeof(tree));
  
    analdef *adef = (analdef*)malloc(sizeof(analdef));   

    /* 
       tell the CPU to ignore exceptions generated by denormalized floating point values.
       If this is not done, depending on the input data, the likelihood functions can exhibit 
       substantial run-time differences for vectors of equal length.
    */
    
#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
    _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
#endif   

  /* get the start time */
   
    masterTime = gettime();         
    
  /* initialize the analysis parameters in struct adef to default values */
    
    initAdef(adef);

  /* parse command line arguments: this has a side effect on tr struct and adef struct variables */
  
    get_args(argc, argv, adef, tr); 
  
  /* generate the ExaML output file names and store them in strings */
    
    makeFileNames();
    
#ifdef _USE_PTHREADS
    /* :TODO: do this differently */
    assert(0); 
    startPthreads(tr);
    
    assert(NumberOfThreads != 0); 
    /* masterBarrier(THREAD_INIT_PARTITION, tr); */
    /* if(!adef->readTaxaOnly) */
    /*   masterBarrier(THREAD_ALLOC_LIKELIHOOD, tr); */
#endif


    initializeTree(tr, adef); 
    
    if(tr->manyPartitions && tr->NumberOfModels < mpiState.commSize && mpiState.rank == 0 )
      {
	printf("You specified -Q for assigning full partitions to processes, but you have less partitions than processes. This would be highly inefficient. Aborting.\n"); 
	MPI_Abort(mpiState.comm, -1);
      }
    
    if(mpiState.rank == 0)  
      {
	printModelAndProgramInfo(tr, adef, argc, argv);  
	printBothOpen("Memory Saving Option: %s\n", (tr->saveMemory == TRUE)?"ENABLED":"DISABLED");   	             
      }  
                         
    /* 
       this will re-start ExaML exactly where it has left off from a checkpoint file,
       while checkpointing is important and has to be implemented for the library we should not worry about this right now 
    */
  
    if(adef->useCheckpoint)
      {      
	/* read checkpoint file */
	restart(tr);       	
	
	/* continue tree search where we left it off */
	computeBIGRAPID(tr, adef, TRUE); 
      }
    else
      {
	/* not important, only used to keep track of total accumulated exec time 
	   when checkpointing and restarts were used */
	
	if(mpiState.rank == 0)
	  accumulatedTime = 0.0;
	
	/* get the starting tree: here we just parse the tree passed via the command line 
	   and do an initial likelihood computation traversal 
	   which we maybe should skeip, TODO */
	       	
	getStartingTree(tr);     
	   	          
	/* 
	   here we do an initial full tree traversal on the starting tree using the Felsenstein pruning algorithm 
	   This should basically be the first call to the library that actually computes something :-)
	*/
      
	evaluateGeneric(tr, tr->start, TRUE);	
	
	/* the treeEvaluate() function repeatedly iterates over the entire tree to optimize branch lengths until convergence */
      	
	treeEvaluate(tr, 1); 

	/* now start the ML search algorithm */
      

	computeBIGRAPID(tr, adef, TRUE); 			     
      }            
      
    /* print some more nonsense into the ExaML_info file */
  
    if(mpiState.rank == 0)
      finalizeInfoFile(tr, adef);
  }
  
  /* return 0 which means that our unix program terminated correctly, the return value is not 1 here */

  MPI_Barrier(mpiState.comm);
  MPI_Finalize();

  return 0;
}


