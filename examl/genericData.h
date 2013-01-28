#ifndef _GENERIC_DATA_H
#define _GENERIC_DATA_H

typedef enum {
  typeUndefined,
  typeDouble, 		
  typeInt, 
} generic_type ; 


typedef struct _G_D_
{
  generic_type aType; 
  
  
  union
  {
    int *intData; 
    double *doubleData; 
  } value; 
    
} GENERIC_DATA ;  


extern   void g_free(GENERIC_DATA *gd); 
extern void g_assignElem(const GENERIC_DATA *dest, int posDest, const GENERIC_DATA *src, int posSrc ) ; 
extern GENERIC_DATA* galloc( int num, generic_type aType); 
extern void g_memcpy(const GENERIC_DATA *dest, int destStart, const GENERIC_DATA* src, int startSrc , int num); 
#endif


