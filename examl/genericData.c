#include "axml.h"


 void g_free(GENERIC_DATA *gd)
{
  switch(gd->aType)
    {
    case typeDouble :
      free(gd->value.doubleData);
      break;
    case typeInt :
      free(gd->value.intData);
      break;
    default: assert(0);
    }

  free(gd);
}



void g_assignElem(const GENERIC_DATA *dest, int posDest, const GENERIC_DATA *src, int posSrc ) 
{
  assert(dest->aType == src->aType) ;
  switch(src->aType)
    {
    case typeDouble:
      dest->value.doubleData[posDest] =  src->value.doubleData[posSrc]; 
      break; 
    case typeInt: 
      dest->value.intData[posDest] =  src->value.intData[posSrc]; 
      break; 
    default: assert(0); 
    }
}


GENERIC_DATA* galloc( int num, generic_type aType)
{
  GENERIC_DATA *result = (GENERIC_DATA*) malloc(sizeof(GENERIC_DATA)); 
  result->aType = aType; 

  switch(aType)
    {
    case typeDouble: 
      result->value.doubleData = (double*)calloc(num, sizeof(double)); 
      break; 
    case typeInt: 
      result->value.intData = (int*)calloc(num, sizeof(int)); 
      break; 
    default :
      assert(0); 
    }

  return result; 
}


void g_memcpy(const GENERIC_DATA *dest, int destStart, const GENERIC_DATA* src, int startSrc , int num)
{  
  switch(dest->aType)
    {
    case typeInt: 
      memcpy(dest->value.intData + destStart,  src->value.intData + startSrc, num * sizeof(int)); 
      break; 
    case typeDouble: 
      memcpy(dest->value.doubleData + destStart, src->value.doubleData + startSrc , num * sizeof(double)); 
      break; 
    default: assert(0); 
    }
}

