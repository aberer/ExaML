#include "axml.h"

#ifndef _TREE_H
#define _TREE_H

boolean setupTree (tree *tr); 
void initializeTree(tree *tr, analdef *adef); 
void initializePartitions(tree *tr, FILE *byteFile); 
int partComparePtr(const void *p1, const void *p2); 


#endif 
