#ifdef _USE_RTS

#ifndef _FAULT_TOLERANCE_H
#define _FAULT_TOLERANCE_H
#include "axml.h" 

#define ACCESS_LIST_INT(elem) *((int*)(elem->value))

void handleMPIError(tree *tr); 

#endif

#else 
void doNothing();
#endif
