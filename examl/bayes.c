#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>



#include "axml.h"
extern char run_id[128];
extern int processID;
extern int Thorough;


typedef enum
{
  SPR,
  stNNI,
  UPDATE_ALL_BL,
  UPDATE_MODEL,
  UPDATE_GAMMA,
}prop;

typedef struct {
  
  /* these 3 are independent of the state, can be taken out unless we want to pass a single pointer as an argument*/  
  nodeptr * list; /* list of possible re-insertion nodes */ 
  int maxradius;  /* maximum radius of re-insertion from the pruning point */
  tree * tr;
  
  double curprior;
  double newprior;
  double hastings;
  
  /*
   * these are the bits that are necessary for the topo moves
   */
  nodeptr p; /* node pruned */
  nodeptr nb;   /* p->next->back describes an edge that dissapears when p is pruned */
  double nbz[NUM_BRANCHES];
  nodeptr nnb; /* p->next->next->back describes an edge that dissapears when p is pruned */
  double nnbz[NUM_BRANCHES];
  nodeptr r; /* edge neighbour of re-insertion node, q->back */
  nodeptr q; /* re-insertion node */
  /*
   * these are the bits that are necessary for the NNI moves
   */
  int whichNNI; /* 0 same topo, 1, 2 */

  /*
   * necessary for the branch length moves
   */
  double qz[NUM_BRANCHES]; /* BL values prior to re-insertion */
  double bl_sliding_window_w;
  double bl_prior;
  double bl_prior_exp_lambda;
  /*
   * necessary for model moves
   */
  analdef * adef;
  int model;
  int nstates;
  int numSubsRates;
  double rt_sliding_window_w;
  double *curSubsRates;//used for resetting

  /*
   * necessary for gamma
   */
  double curAlpha;
  double gm_sliding_window_w;

} state;


static state *state_init(tree *tr, analdef * adef, int maxradius, double bl_w, double rt_w, double gm_w, double bl_p)
{
  state *curstate  =(state *)malloc(sizeof(state));
  nodeptr *list = (nodeptr *)malloc(sizeof(nodeptr) * 2 * tr->mxtips);
  curstate->list = list;
  curstate->maxradius = maxradius;
  curstate->tr = tr;
  curstate->bl_sliding_window_w = bl_w;
  curstate->bl_prior = 1.0;
  curstate->bl_prior_exp_lambda = bl_p;
  //this can be extended to more than one partition, but one for now
  curstate->model = 0;
  curstate->adef = adef;
  curstate->rt_sliding_window_w = rt_w;
  curstate->nstates = tr->partitionData[curstate->model].states; /* 4 for DNA */
  curstate->numSubsRates = (curstate->nstates * curstate->nstates - curstate->nstates) / 2; /* 6 for DNA */
  curstate->curSubsRates = (double *) malloc(curstate->numSubsRates * sizeof(double));
  curstate->gm_sliding_window_w = gm_w;
  assert(curstate != NULL);
  return curstate;
}


static void printSubsRates(tree *tr,int model, int numSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  printBothOpen("Subs rates: ");
  for(i=0; i<numSubsRates; i++)
    printBothOpen("%d => %.3f, ", i, tr->partitionData[model].substRates[i]);
  printBothOpen("\n\n");
}

static void recordSubsRates(tree *tr, int model, int numSubsRates, double *prevSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    prevSubsRates[i] = tr->partitionData[model].substRates[i];
}


static void reset_branch_length(nodeptr p, int numBranches)
{
  int i;
  double new_value;
  for(i = 0; i < numBranches; i++)
  {
    assert(p->z_tmp[i] == p->back->z_tmp[i]);
    p->z[i] = p->back->z[i] = p->z_tmp[i];   /* restore saved value */
  }
}

static double exp_pdf(double lambda, double x)
{
  return (lambda * exp(-(lambda * x))); 
}

//setting this out to allow for other types of setting
static void set_branch_length_sliding_window(nodeptr p, int numBranches,state * s, boolean record_tmp_bl)
{
  int i;
  double new_value;
  double r,mx,mn;
  for(i = 0; i < numBranches; i++)
  {
    double real_z;
    
    if(record_tmp_bl)
    {
      assert(p->z[i] == p->back->z[i]); 
      p->z_tmp[i] = p->back->z_tmp[i] = p->z[i];   /* keep current value */
    }
    r = (double)rand()/(double)RAND_MAX;
    
    
    
    real_z = -log(p->z[i]) * s->tr->fracchange;
    
//     printf( "z: %f %f\n", p->z[i], real_z );
    
    mn = real_z-(s->bl_sliding_window_w/2);
    mx = real_z+(s->bl_sliding_window_w/2);
    new_value = exp(-(fabs(mn + r * (mx-mn)/s->tr->fracchange )));
    
    
    /* Ensure always you stay within this range */
    if(new_value > zmax) new_value = zmax;
    if(new_value < zmin) new_value = zmin;
    assert(new_value <= zmax && new_value >= zmin);
//     printf( "z: %f %f %f\n", p->z[i], new_value, real_z );
    p->z[i] = p->back->z[i] = new_value;
    //assuming this will be visiting each node, and multiple threads won't be accessing this
    //s->bl_prior += log(exp_pdf(s->bl_prior_exp_lambda,new_value));
    //s->bl_prior += 1;
  }
}

static void traverse_branches(nodeptr p, int *count, state * s, boolean resetBL)
{
  nodeptr q;
  //printf("current BL at %db%d: %f\n", p->number, p->back->number, p->z[0]);
  if(resetBL)
    reset_branch_length(p, s->tr->numBranches);
  else//can allow for other methods later
    set_branch_length_sliding_window(p, s->tr->numBranches, s, TRUE);
  *count += 1;


  if (! isTip(p->number, s->tr->mxtips)) 
  {                                  /*  Adjust descendants */
    q = p->next;
    while (q != p) 
    {
      traverse_branches(q->back, count, s, resetBL);
      q = q->next;
    }   
    newviewGeneric(s->tr, p, FALSE);     // not sure if we need this
  }
}


static void update_all_branches(state * s, boolean resetBL)
{
  int updated_branches = 0;
  assert(isTip(s->tr->start->number, s->tr->mxtips));
  /* visit each branch exactly once */
  traverse_branches(s->tr->start->back, &updated_branches, s, resetBL);
  assert(updated_branches == s->tr->mxtips + s->tr->mxtips - 3);
}


static void traverse_branches_set_fixed(nodeptr p, int *count, state * s, double z )
{
  nodeptr q;
  int i;
  //printf("current BL at %db%d: %f\n", p->number, p->back->number, p->z[0]);
  
  for( i = 0; i < s->tr->numBranches; i++)
  {
   
    p->z[i] = p->back->z[i] = z;
  }
  
  *count += 1;


  if (! isTip(p->number, s->tr->mxtips)) 
  {                                  /*  Adjust descendants */
    q = p->next;
    while (q != p) 
    {
      traverse_branches_set_fixed(q->back, count, s, z);
      q = q->next;
    }   
    newviewGeneric(s->tr, p, FALSE);     // not sure if we need this
  }
}


/*
 * should be sliding window proposal
 */

static boolean simpleBranchLengthProposal(state * instate)
{
   
  //for each branch get the current branch length
  //pull a uniform like
  //x = current, w =window
  //uniform(x-w/2,x+w/2)

  update_all_branches(instate, FALSE);
  evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* update the tr->likelihood */

  //for prior, just using exponential for now
  //calculate for each branch length
  // where lambda is chosen and x is the branch length
  //lambda * exp(-lamba * x)

  //only calculate the new ones
  //
  return TRUE;
}

static void resetSimpleBranchLengthProposal(state * instate)
{
  update_all_branches(instate, TRUE);
}


/*
 * should be sliding window proposal
 */

static void editSubsRates(tree *tr, int model, int subRatePos, double subRateValue)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  assert(subRateValue <= RATE_MAX && subRateValue >= RATE_MIN);
  int states = tr->partitionData[model].states; 
  int numSubsRates = (states * states - states) / 2;
  assert(subRatePos >= 0 && subRatePos < numSubsRates);
  tr->partitionData[model].substRates[subRatePos] = subRateValue;
}

static void simpleModelProposal(state * instate)
{
  //TODO: add safety to max and min values
  //record the old ones
  recordSubsRates(instate->tr, instate->model, instate->numSubsRates, instate->curSubsRates);
  //choose a random set of model params,
  //probably with dirichlet proposal
  //with uniform probabilities, no need to have other
  int state;
  double new_value,curv;
  double r,mx,mn;
  //using the branch length sliding window for a test
  for(state = 0;state<instate->numSubsRates ; state ++)
    {
      curv = instate->tr->partitionData[instate->model].substRates[state];
      r = (double)rand()/(double)RAND_MAX;
      mn = curv-(instate->rt_sliding_window_w/2);
      mx = curv+(instate->rt_sliding_window_w/2);
      new_value = fabs(mn + r * (mx-mn));
      /* Ensure always you stay within this range */
      if(new_value > RATE_MAX) new_value = RATE_MAX;
      if(new_value < RATE_MIN) new_value = RATE_MIN;
      //printf("%i %f %f\n", state, curv, new_value);
      editSubsRates(instate->tr,instate->model, state, new_value);
    }
  //recalculate eigens

  initReversibleGTR(instate->tr, instate->model); /* 1. recomputes Eigenvectors, Eigenvalues etc. for Q decomp. */

  /* TODO: need to broadcast rates here for parallel version ! */

  evaluateGeneric(instate->tr, instate->tr->start, TRUE); /* 2. re-traverse the full tree to update all vectors */
  //TODO: without this, the run will fail after a successful model, but failing SPR
  //TODOFER: what did we have in mind regarding the comment above?
  
  evaluateGeneric(instate->tr, instate->tr->start, FALSE);
  //for prior, just use dirichlet
  // independent gamma distribution for each parameter
  //the pdf for this is
  // for gamma the prior is gamma

  //for statefreqs should all be uniform

  //only calculate the new ones
}

static void restoreSubsRates(tree *tr, analdef *adef, int model, int numSubsRates, double *prevSubsRates)
{
  assert(tr->partitionData[model].dataType = DNA_DATA);
  int i;
  for(i=0; i<numSubsRates; i++)
    tr->partitionData[model].substRates[i] = prevSubsRates[i];

  initReversibleGTR(tr, model);

  /* TODO need to broadcast rates here for parallel version */

  evaluateGeneric(tr, tr->start, TRUE);
}

static void resetSimpleModelProposal(state * instate)
{
  restoreSubsRates(instate->tr, instate->adef, instate->model, instate->numSubsRates, instate->curSubsRates);
  //evaluateGeneric(instate->tr, instate->tr->start, FALSE);
}

static nodeptr selectRandomSubtree(tree *tr)
{
  nodeptr 
    p;

  do
    {
      int 
        exitDirection = rand() % 3; 
     
      p = tr->nodep[(rand() % (tr->mxtips - 2)) + 1 + tr->mxtips];
      
      switch(exitDirection)
        {
        case 0:
          break;
        case 1:
          p = p->next;
          break;
        case 2:
          p = p->next->next;
          break;
        default:
          assert(0);
        }
    }
  while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));

  assert(!isTip(p->number, tr->mxtips));

  return p;
}

static void recordBranchInfo(nodeptr p, double *bl, int numBranches)
{
  int i;
  for(i = 0; i < numBranches; i++)
    bl[i] = p->z[i];
}
static int spr_depth = 0;

static node *randomSPR_traverse( tree *tr, node *n ) {

  double randprop = (double)rand()/(double)RAND_MAX;
  
  if( isTip(n->number, tr->mxtips ) || randprop < 0.1 ) {
    return n;
  } else if( randprop >= 0.1 && randprop < 0.55 ) {
    spr_depth++;
    return randomSPR_traverse( tr, n->next->back );
  } else {
    spr_depth++;
    return randomSPR_traverse( tr, n->next->next->back );
  }
}
static node *randomSPR( tree *tr, node *n ) {
  if( isTip(n->number, tr->mxtips ) ) {
    return n;
  }
  spr_depth = 1;
  double randprop = (double)rand()/(double)RAND_MAX;
  if( randprop < 0.5 ) {
    return randomSPR_traverse( tr, n->next->back );
  } else {
    return randomSPR_traverse( tr, n->next->next->back );
  }
}

static void doSPR(tree *tr, state *instate)
{
  nodeptr    
    p = selectRandomSubtree(tr);
  
  /* evaluateGeneric(tr, tr->start, TRUE);
     printf("%f \n", tr->likelihood);*/

#if 0
  parsimonySPR(p, tr);
#endif

  

  /*evaluateGeneric(tr, tr->start, TRUE);
    printf("%f \n", tr->likelihood);*/

  instate->p = p;
  instate->nb  = p->next->back;
  instate->nnb = p->next->next->back;
  
  recordBranchInfo(instate->nb, instate->nbz, instate->tr->numBranches);
  recordBranchInfo(instate->nnb, instate->nnbz, instate->tr->numBranches);

  /* removeNodeBIG(tr, p,  tr->numBranches); */
  /* remove node p */
  double   zqr[NUM_BRANCHES];
  int i;
  for(i = 0; i < tr->numBranches; i++)
  {
    zqr[i] = instate->nb->z[i] * instate->nnb->z[i];        
    if(zqr[i] > zmax) zqr[i] = zmax;
    if(zqr[i] < zmin) zqr[i] = zmin;
  }
  hookup(instate->nb, instate->nnb, zqr, tr->numBranches); 
  p->next->next->back = p->next->back = (node *) NULL;
  /* done remove node p (omitted BL opt) */

  
  double randprop = (double)rand()/(double)RAND_MAX;
  spr_depth = 0;
  assert( !(isTip(instate->nb->number, instate->tr->mxtips ) && isTip(instate->nnb->number, instate->tr->mxtips )) );
  if( isTip(instate->nb->number, instate->tr->mxtips ) ) {
    randprop = 1.0;
  } else if( isTip(instate->nnb->number, instate->tr->mxtips ) ) {
    randprop = 0.0;
  }
  
  if( randprop <= 0.5 ) {
    tr->insertNode = randomSPR( tr, instate->nb );
  } else {
    tr->insertNode = randomSPR( tr, instate->nnb );
  }

  

  instate->q = tr->insertNode;
  instate->r = instate->q->back;
  recordBranchInfo(instate->q, instate->qz, instate->tr->numBranches);

  
  Thorough = 0;
 // assert(tr->thoroughInsertion == 0);
  /* insertBIG wont change the BL if we are not in thorough mode */
  
  insertBIG(instate->tr, instate->p, instate->q, instate->tr->numBranches);
  evaluateGeneric(instate->tr, instate->p->next->next, FALSE); 
  /*testInsertBIG(tr, p, tr->insertNode);*/

  //printf("%f \n", tr->likelihood);
}

static void resetSPR(state * instate)
{
  /* prune the insertion */
  hookup(instate->q, instate->r, instate->qz, instate->tr->numBranches);
  instate->p->next->next->back = instate->p->next->back = (nodeptr) NULL;
  /* insert the pruned tree in its original node */
  hookup(instate->p->next,       instate->nb, instate->nbz, instate->tr->numBranches);
  hookup(instate->p->next->next, instate->nnb, instate->nnbz, instate->tr->numBranches);
  newviewGeneric(instate->tr, instate->p, FALSE); 
}


static prop proposal(state * instate)
/* so here the idea would be to randomly choose among proposals? we can use typedef enum to label each, and return that */ 
{
  double randprop = (double)rand()/(double)RAND_MAX;
  boolean proposalSuccess;
  //double start_LH = evaluateGeneric(instate->tr, instate->tr->start); /* for validation */
  prop proposal_type;
  //simple proposal
  
#if 0
  if(randprop < 0.25) 
  {
    if(randprop < 0.2)//TOPOLOGICAL MOVE
    {
      if(randprop > 0.1)//SPR MOVE
      {
        proposal_type = SPR;
        if (randprop < 0.15)
          instate->maxradius = 1;
        else
          instate->maxradius = 2;

        doSPR(instate->tr, instate);

        proposalSuccess = TRUE;

        /* TODO  shall we still want to use this as an alternative to the parsimony-informed SPR?*/
        /*proposalSuccess = simpleNodeProposal(instate);*/
      }
      else
      {
        proposal_type = stNNI;
        proposalSuccess = stNNIproposal(instate); 
        if(proposalSuccess == FALSE)
        {
          printBothOpen("ERROR: stNNI proposal failed\n");
          assert(FALSE);
        }
      }
      if(proposalSuccess == FALSE)
      {
        assert(FALSE); // this should either never happen or look below and return PROPOSAL_FAILED to react accordingly
      }
      else
      {
        /* A moved has been made, previous state is in instate */
        if(proposal_type != stNNI) /*TODOFER delete this when bl are changed (bl should change always in the stNNI?)*/
          assert(instate->tr->startLH != instate->tr->likelihood);
      }
    }
    else{//MODEL
      proposal_type = UPDATE_MODEL;
      simpleModelProposal(instate);
    }
  }
  else
  {
    if(randprop < 0.95)//UPDATE_ALL_BL
    {
      proposal_type = UPDATE_ALL_BL;
      instate->bl_prior = 0;
      //printBothOpen("Propose BL_UPDATE\n");
      assert(proposal_type == UPDATE_ALL_BL);
      proposalSuccess = simpleBranchLengthProposal(instate);
      assert(instate->tr->startLH != instate->tr->likelihood);
      assert(proposalSuccess);
    }
    else//GAMMA
    {
      proposal_type = UPDATE_GAMMA;
      simpleGammaProposal(instate);
    }
  }
  //record the curprior
  instate->newprior = instate->bl_prior;
  //print_proposal(proposal_type);
  return proposal_type;
#else
//   printf( "%d %f\n", processID, randprop );
  instate->newprior = instate->bl_prior;
  if( randprop < 0.1 ) {
    proposal_type = UPDATE_MODEL;
    simpleModelProposal(instate);
  } else if( 1 && randprop < 0.2 ) {
    proposal_type = UPDATE_ALL_BL;
    proposalSuccess = simpleBranchLengthProposal(instate);
    
  } else {
    
    proposal_type = SPR;
    
    if (randprop < 0.15)
      instate->maxradius = 1;
    else
      instate->maxradius = 2;
    
    doSPR(instate->tr, instate);
    
//      printf( "propose SPR: %d\n", spr_depth );
  }
  return proposal_type;
#endif
  
}

static void resetState(prop proposal_type, state * curstate)
{
  switch(proposal_type)
  {

    case SPR:
      resetSPR(curstate);
      break;
#if 0
    case stNNI:
      reset_stNNI(curstate);
      break;
#endif
    case UPDATE_ALL_BL:
      resetSimpleBranchLengthProposal(curstate);
      //printf("RESETBL\n");
      break;

    case UPDATE_MODEL:
      resetSimpleModelProposal(curstate);
      break;
#if 0
    case UPDATE_GAMMA:
      resetSimpleGammaProposal(curstate);
      break;
#endif
    default:
      assert(FALSE);
  }
}

#if 0
static void printRecomTree(tree *tr, boolean printBranchLengths, char *title)
{
  FILE *nwfile;
  nwfile = myfopen("tmp.nw", "w+");
  Tree2StringRecomREC(tr->tree_string, tr, tr->start->back, printBranchLengths);
  fprintf(nwfile,"%s\n", tr->tree_string);
  fclose(nwfile);
  if(title)
    printBothOpen("%s\n", title);
  if (printBranchLengths)
    printBothOpen("%s\n", tr->tree_string);
  printBothOpen("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  //system("bin/nw_display tmp.nw");
}  
#endif


static void printStateFile(int iter, state * curstate)
{ 
  const size_t tmp_len = 256;
  char tmp[tmp_len];
  
  strncpy(tmp, "RAxML_states.", tmp_len);
  strncat(tmp, run_id, tmp_len);
  FILE *f = myfopen(tmp, "ab");
  fprintf(f,"%d\t%f",iter, curstate->tr->likelihood);
  int i;
  for(i = 0;i < curstate->numSubsRates; i++)
  {
    fprintf(f,"\t%f",curstate->curSubsRates[i]);
  }
  fprintf(f,"\t%f",curstate->curAlpha);
  fprintf(f,"\n");
  fclose(f);
}


node *find_tip( node *n, tree *tr ) {
  if( isTip(n->number, tr->mxtips) ) {
    return n;
  } else {
    return find_tip( n->back, tr );
  }
  
}
static char *Tree2StringRecomREC(char *treestr, tree *tr, nodeptr q, boolean printBranchLengths)
{
  char  *nameptr;            
  double z;
  nodeptr p = q;
  int slot;

  if(isTip(p->number, tr->mxtips)) 
  {               
    nameptr = tr->nameList[p->number];     
    sprintf(treestr, "%s", nameptr);
    while (*treestr) treestr++;
  }
  else 
  {                      
    while(!p->x)
      p = p->next;
    *treestr++ = '(';
    treestr = Tree2StringRecomREC(treestr, tr, q->next->back, printBranchLengths);
    *treestr++ = ',';
    treestr = Tree2StringRecomREC(treestr, tr, q->next->next->back, printBranchLengths);
    if(q == tr->start->back) 
    {
      *treestr++ = ',';
      treestr = Tree2StringRecomREC(treestr, tr, q->back, printBranchLengths);
    }
    *treestr++ = ')';                    
    // write innernode as nodenum_b_nodenumback
#if 0
  sprintf(treestr, "%d", q->number);
    while (*treestr) treestr++;
    *treestr++ = 'b';                    
    sprintf(treestr, "%d", p->back->number);
    while (*treestr) treestr++;
#endif
    
  }

  if(q == tr->start->back) 
  {              
    if(printBranchLengths)
      sprintf(treestr, ":0.0;\n");
    else
      sprintf(treestr, ";\n");                  
  }
  else 
  {                   
    if(printBranchLengths)          
    {
      //sprintf(treestr, ":%8.20f", getBranchLength(tr, SUMMARIZE_LH, p));                 
      assert(tr->fracchange != -1.0);
      z = q->z[0];
      if (z < zmin) 
        z = zmin;        
      sprintf(treestr, ":%8.20f", -log(z) * tr->fracchange);               
    }
    else            
      sprintf(treestr, "%s", "\0");         
  }

  while (*treestr) treestr++;
  return  treestr;
}

static int tree_dump_num = 0;
static void printRecomTree(tree *tr, boolean printBranchLengths, char *title)
{
  FILE *nwfile;
  const size_t tmp_len = 256;
  char tmp[tmp_len];
  
  snprintf( tmp, tmp_len, "spr_%04d.txt", tree_dump_num );
  ++tree_dump_num;
  
  nwfile = myfopen(tmp, "w");
  Tree2StringRecomREC(tr->tree_string, tr, tr->start->back, printBranchLengths);
  fprintf(nwfile,"%s\n", tr->tree_string);
  fclose(nwfile);
  if(title)
    printBothOpen("%s\n", title);
  if (printBranchLengths)
    printBothOpen("%s\n", tr->tree_string);
  printBothOpen("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  //system("bin/nw_display tmp.nw");
}  

void mcmc(tree *tr, analdef *adef)
{
  boolean proposalAccepted;
  
  int j;
  int maxradius = 30;
  int num_moves = 100000;
  
  size_t accepted_nni = 0;
  size_t accepted_spr = 0;
  size_t accepted_bl = 0;
  size_t accepted_gamma = 0;
  size_t accepted_model = 0;

  size_t rejected_nni = 0;
  size_t rejected_spr = 0;
  size_t rejected_bl = 0;
  size_t rejected_gamma = 0;
  size_t rejected_model = 0;

  size_t inserts = 0;
  double printTime = 0;
  
   //allocate states
  double bl_prior_exp_lambda = 0.1;
  double bl_sliding_window_w = 0.005;
  double gm_sliding_window_w = 0.75;
  double rt_sliding_window_w = 0.5;
  
  srand(123);
  
  tr->start = find_tip(tr->start, tr );
  
  printf( "isTip: %d %d\n", tr->start->number, tr->mxtips );
  printf( "z: %f\n", tr->start->z[0] );
  
  assert( isTip(tr->start->number, tr->mxtips ));
  
  state *curstate = state_init(tr, adef, maxradius, bl_sliding_window_w, rt_sliding_window_w, gm_sliding_window_w, bl_prior_exp_lambda);
  
  evaluateGeneric(tr, tr->start, TRUE);  
  tr->startLH = tr->likelihood;
  printBothOpen("Starting with tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);

  /* Set reasonable model parameters */
  evaluateGeneric(curstate->tr, curstate->tr->start, FALSE); // just for validation 
  printBothOpen("tr LH before modOpt %f\n",curstate->tr->likelihood);
  printSubsRates(curstate->tr, curstate->model, curstate->numSubsRates);

#if 0
  /* optimize the model with Brents method for reasonable starting points */
  modOpt(curstate->tr, 5.0); /* not by proposal, just using std raxml machinery... */
  evaluateGeneric(curstate->tr, curstate->tr->start, FALSE); // just for validation 
  printBothOpen("tr LH after modOpt %f\n",curstate->tr->likelihood);
  printSubsRates(curstate->tr, curstate->model, curstate->numSubsRates);
  recordSubsRates(curstate->tr, curstate->model, curstate->numSubsRates, curstate->curSubsRates);
#endif

  //printf( "zx: %f\n", tr->start->z[0] );
  evaluateGeneric(tr, tr->start, TRUE);
  printf( "at start: %f\n", tr->likelihood );
  
  if( 0 )
  {
    
    int count = 0;
    traverse_branches_set_fixed( tr->start, &count, curstate, 0.65 );
  }
  
  evaluateGeneric(tr, tr->start, TRUE);
  printf( "after reset start: %f\n", tr->likelihood );
  int first = 1;
  
  curstate->curprior = 1;
  curstate->hastings = 1;
  
  /* beginning of the MCMC chain */
  for(j=0; j<num_moves; j++)
  {
    prop which_proposal;
    double t = gettime(); 
    double proposalTime = 0.0;
    double testr;
    double acceptance;
    
//     printBothOpen("iter %d, tr LH %f, startLH %f\n",j, tr->likelihood, tr->startLH);
    proposalAccepted = FALSE;
    

    evaluateGeneric(tr, tr->start, FALSE); // just for validation (make sure we compare the same)
//      printBothOpen("before proposal, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
    tr->startLH = tr->likelihood;

    which_proposal = proposal(curstate);
//     printf( "z: %f\n", tr->start->z[0] );
    if (first == 1)
    {
      first = 0;
      curstate->curprior = curstate->newprior;
    }
//     printBothOpen("proposal done, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
    assert(which_proposal == SPR || which_proposal == stNNI ||
           which_proposal == UPDATE_ALL_BL || 
           which_proposal == UPDATE_MODEL || which_proposal == UPDATE_GAMMA);
    proposalTime += gettime() - t;
    /* decide upon acceptance */
    testr = (double)rand()/(double)RAND_MAX;
    //should look something like 
    acceptance = fmin(1,(curstate->hastings) * 
                       (exp(curstate->newprior-curstate->curprior)) * (exp(curstate->tr->likelihood-curstate->tr->startLH)));
    
//     acceptance = exp(curstate->tr->likelihood-curstate->tr->startLH);
    
//     printf( "exp: %f\n", (curstate->tr->likelihood-curstate->tr->startLH) );
//     printf( "curstate: %f %f %f\n", curstate->hastings, curstate->newprior, curstate->curprior );
//     printf( "acc: %f %f\n", testr, acceptance );
    /*
      //printRecomTree(tr, FALSE, "after proposal");
      printBothOpen("after proposal, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
    */
    
    
    if(processID == 0 && (j % 100) == 0) {
      printf( "propb: %d %f %f %d %d %d %d %f %f %f\n", j, tr->likelihood, tr->startLH, testr < acceptance, accepted_spr, accepted_model, accepted_bl, curstate->hastings, curstate->newprior, curstate->curprior );
      
      printSubsRates(curstate->tr, curstate->model, curstate->numSubsRates);
    }
    
    if(testr < acceptance)
    {
      proposalAccepted = TRUE;

      switch(which_proposal)
        {
          

        case SPR:      
//           printf( "accept spr: %d\n", spr_depth );
          
#if 0
          if( processID == 0 ) {
            printRecomTree(tr, TRUE, "after accepted");
          }
#endif
          // printBothOpen("SPR new topology , iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
          accepted_spr++;
          break;
        case stNNI:       
          //printBothOpen("NNI new topology , iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
          accepted_nni++;
          break;

        case UPDATE_ALL_BL:       
          //      printBothOpen("BL new , iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
          accepted_bl++;
          break;
        case UPDATE_MODEL:      
          //    printBothOpen("Model new, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
          accepted_model++;
          break;
        case UPDATE_GAMMA:      
          //    printBothOpen("Gamma new, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
          accepted_gamma++;
          break;
        default:
          assert(0);
        }

      //printBothOpen("accepted , iter %d tr LH %f, startLH %f, %i \n", j, tr->likelihood, tr->startLH, which_proposal);
      curstate->tr->startLH = curstate->tr->likelihood;  //new LH
      curstate->curprior = curstate->newprior;          
    }
    else
    {
      //printBothOpen("rejected , iter %d tr LH %f, startLH %f, %i \n", j, tr->likelihood, tr->startLH, which_proposal);
    //print_proposal(which_proposal);
      resetState(which_proposal,curstate);
      
      switch(which_proposal)
        {
        case SPR:
          rejected_spr++;
          break;
        case stNNI:
          rejected_nni++;
          break;
        case UPDATE_ALL_BL:
          rejected_bl++;
          break;
        case UPDATE_MODEL:
          rejected_model++;
          break;
        case UPDATE_GAMMA:
          rejected_gamma++;
          break;
        default:
          assert(0);
        }
      
      evaluateGeneric(tr, tr->start, FALSE); 
      
//       printf( "accept: %d\n", proposalAccepted );
      
      
      // just for validation 

      if(fabs(curstate->tr->startLH - tr->likelihood) > 1.0E-15)
      {
//         printRecomTree(tr, TRUE, "after reset");
        printBothOpen("WARNING: LH diff %.20f\n", curstate->tr->startLH - tr->likelihood);
        printBothOpen("after reset, iter %d tr LH %f, startLH %f\n", j, tr->likelihood, tr->startLH);
      }
      //assert(fabs(curstate->tr->startLH - tr->likelihood) < 1.0E-10); 
      assert(fabs(curstate->tr->startLH - tr->likelihood) < 0.1);
    }       
    inserts++;
    
#if 0
    /* need to print status */
    if (j % 50 == 0)
    {
      t = gettime(); 
      printBothOpen("sampled at iter %d, tr LH %f, startLH %f, prior %f, incr %f\n",j, tr->likelihood, tr->startLH, curstate->curprior, tr->likelihood - tr->startLH);
      boolean printBranchLengths = TRUE;
      /*printSimpleTree(tr, printBranchLengths, adef);*/
      //TODO: print some parameters to a file 
      printStateFile(j,curstate);
      printTime += gettime() - t;
    }
#endif
    
  }

}