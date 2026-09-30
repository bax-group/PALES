/******************************************************************************/
/*                   ---- NIH NMR Software System ----                        */
/*                        Copyright 1992 to 2015                              */
/*                             Frank Delaglio                                 */
/*                   NIH Laboratory of Chemical Physics                       */
/*                                                                            */
/*               This software is not for distribution without                */
/*                  the written permission of the author.                     */
/******************************************************************************/

/***/
/* memory.c
/*
/* Procedure wrappers to allocate dynamic memory, for
/* monitoring and dubugging of memory use.
/*
/* Many of these procedures have wrappers in memory.h
/* to recast argumments to proper types.
/***/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "memory.h"

static  NMR_INT bytesUsed = 0, bytesReused = 0;
static  int verboseAlloc = 0;

#define FPR (void)fprintf

#define ERRTEST( P ) \
   if (!(P)) (void)fprintf( stderr, "Memory allocation failure.\n" );

/***/
/* Adjust verbose flag:
/***/

int setAlloc( verboseFlag )

   int verboseFlag;
{
   return( verboseAlloc = verboseFlag );
}

/***/
/* Return memory in use:
/***/

int getMemUsed( usedPtr, reusedPtr )

   NMR_INT *usedPtr, *reusedPtr;
{

   *usedPtr   = bytesUsed;
   *reusedPtr = bytesReused;

   return( 0 );
}

/***/
/* Show current memory usage.
/***/

int showAlloc()
{

#ifdef NMR64
    FPR( stderr, "REMARK Memory Allocation Status. Reused: %ld In Use: %ld\n", bytesReused, bytesUsed );
#else
    FPR( stderr, "REMARK Memory Allocation Status. Reused: %d In Use: %d\n", bytesReused, bytesUsed );
#endif

    return( 0 );
}

/***/
/* Show current memory use.
/***/

int ShowCurrentAlloc64( type, caller, ptr, length )

   char *type, *caller;
   void *ptr;
   NMR_INT length;
{
#ifdef NMR64
   FPR( stderr, "\n_MEM %s Name: %s Alloc: %p Add: %ld Now: %ld\n", type, caller, ptr, length, bytesUsed );
#else
   FPR( stderr, "\n_MEM %s Name: %s Alloc: %p Add: %d Now: %d\n", type, caller, ptr, length, bytesUsed );
#endif
   return( 0 );
}

/***/
/* Allocate a float array:
/***/

float *fltMalloc64( caller, length )

   char    *caller;
   NMR_INT length;
{
   float *ptr;

   length    *= sizeof(float);
   bytesUsed += length;
   ptr        = (float *) malloc( (size_t)length );

   if (verboseAlloc)
      {
       (void) showCurrentAlloc( "FLT_ALLOC", caller, ptr, length );
      }

   ERRTEST( ptr );
 
   return( ptr );
}

/***/
/* Allocate an int array:
/***/

int *intMalloc64( caller, length )

   char    *caller;
   NMR_INT length;
{
   int *ptr;

   length    *= sizeof(int);
   bytesUsed += length;
   ptr        = (int *) malloc( (size_t)length );

   if (verboseAlloc)
      {
       (void) showCurrentAlloc( "INT_ALLOC", caller, ptr, length );
      }

   ERRTEST( ptr );

   return( ptr );
}

/***/
/* Allocate an NMR_INT array:
/***/

NMR_INT *n64Malloc( caller, length )

   char    *caller;
   NMR_INT length;
{
   NMR_INT *ptr;

   length    *= sizeof(NMR_INT);
   bytesUsed += length;
   ptr        = (NMR_INT *) malloc( (size_t)length );

   if (verboseAlloc)
      {
       (void) showCurrentAlloc( "NMRINT_ALLOC", caller, ptr, length );
      }

   ERRTEST( ptr );

   return( ptr );
}

/***/
/* Allocate a double prec float array:
/***/

double *dFltMalloc64( caller, length )

   char    *caller;
   NMR_INT length;
{
   double *ptr;

   ptr        = (double *)NULL;
   length    *= sizeof(double);
   bytesUsed += length;

   ptr        = (double *) malloc( (size_t)length );

   if (verboseAlloc)
      {
       (void) showCurrentAlloc( "DBL_ALLOC", caller, ptr, length );
      }

   ERRTEST( ptr );

   return( ptr );
}

/***/
/* Allocate a long int array:
/***/

long int *dIntMalloc64( caller, length )

   char    *caller;
   NMR_INT length;
{
   long int *ptr;

   length    *= sizeof(long int);
   bytesUsed += length;
   ptr        = (long int *) malloc( (size_t)length );

   if (verboseAlloc)
      {
       (void) showCurrentAlloc( "LNG_ALLOC", caller, ptr, length );
      }

   ERRTEST( ptr );

   return( ptr );
}

/***/
/* Allocate a character array:
/***/

char *charMalloc64( caller, length )

   char    *caller;
   NMR_INT length;
{
   char *ptr;

   length    *= sizeof(char);
   bytesUsed += length;
   ptr        = (char *) malloc( (size_t)length );

   if (verboseAlloc)
      {
       (void) showCurrentAlloc( "CHR_ALLOC", caller, ptr, length );
      }

   ERRTEST( ptr );

   return( ptr );
}

/***/
/* strFree: deallocate string created by strDup.
/***/

int strFree( string )

   char *string;
{
    if (!string) return( 0 );

    (void) deAlloc( "strDup", string, 1+strlen(string) );

    return( 0 );
}

/***/
/* strDup: allocate, create, and return a duplicate of string.
/***/

char *strDup( string )

   char *string;
{
    char *dup;

    if (!(dup = charAlloc( "strDup", 1+strlen(string) )))
       {
        ERRTEST( dup );
        return( dup );
       }

    (void) strcpy( dup, string );

    return( dup );
}

/***/
/* strDup2: allocate, create, and return a duplicate of string.
/***/

char *strDup2( caller, string )

   char *caller, *string;
{
    char *dup;

    if (!(dup = charAlloc( caller, 1+strlen(string) )))
       {
        ERRTEST( dup );
        return( dup );
       }

    (void) strcpy( dup, string );

    return( dup );
}

/***/
/* Allocate a type void array:
/***/

void *voidMalloc64( caller, length )

   char    *caller;
   NMR_INT length;
{
   void *ptr;

   bytesUsed += length;
   ptr        = (void *) malloc( (size_t)length );

   if (verboseAlloc)
      {
       (void) showCurrentAlloc( "VOI_ALLOC", caller, ptr, length );
      }

   ERRTEST( ptr );

   return( ptr );
}

/***/
/* Allocate a float matrix.
/***/

float **MatMalloc64( caller, xLength, yLength )

   char    *caller;
   NMR_INT xLength, yLength;
{
    float *work, **matrix;

    bytesUsed += sizeof(float)*xLength*yLength + sizeof(float *)*yLength;

    if (!(work = fltAlloc( caller, yLength*xLength )))
       {
        return( (float **) NULL );
       }

    if (!(matrix = Mat2Ptr64( work, xLength, yLength )))
       {
        (void) deAlloc( caller, work, sizeof(float)*yLength*xLength );
        return( (float **) NULL );
       }

    return( matrix );
}

/***/
/* Allocate a double matrix.
/***/

double **MatMallocD64( caller, xLength, yLength )

   char    *caller;
   NMR_INT xLength, yLength;
{
    double *work, **matrix;

    bytesUsed += sizeof(float)*xLength*yLength + sizeof(float *)*yLength;

    if (!(work = dFltAlloc( caller, yLength*xLength )))
       {
        return( (double **) NULL );
       }

    if (!(matrix = Mat2PtrD64( work, xLength, yLength )))
       {
        (void) deAlloc( caller, work, sizeof(float)* yLength*xLength );
        return( (double **) NULL );
       }

    return( matrix );
}

/***/
/* Reallocate a type int array:
/***/

int *intReMalloc64( caller, ptr, newLength, oldLength )

   NMR_INT newLength, oldLength;
   char    *caller;
   int     *ptr;
{
   int *newPtr;

   if (!ptr) oldLength = 0;

   oldLength   *= sizeof(int);
   newLength   *= sizeof(int);

   bytesUsed   += newLength - oldLength;
   bytesReused += oldLength;

   if (verboseAlloc)
      {
       FPR( stderr,
            "\n_MEM INT_RE_AL Name: %s Free: %p",
            caller, ptr );
      }

   if (ptr)
      newPtr = (int *) realloc( (void *)ptr, (size_t)newLength );
   else
      newPtr = (int *) malloc( (size_t)newLength );

   if (verboseAlloc)
      {
       FPR( stderr, " Alloc: %p", newPtr );

#ifdef NMR64
       FPR( stderr, " Add: %ld Delete: %ld Now: %ld\n", newLength, oldLength, bytesUsed );
#else
       FPR( stderr, " Add: %d Delete: %d Now: %d\n", newLength, oldLength, bytesUsed );
#endif
      }

   ERRTEST( newPtr );

   return( newPtr );
}

/***/
/* Reallocate a type float array:
/***/

float *fltReMalloc64( caller, ptr, newLength, oldLength )

   NMR_INT newLength, oldLength;
   char    *caller;
   float   *ptr;
{
   float *newPtr;

   if (!ptr) oldLength = 0;

   oldLength   *= sizeof(float);
   newLength   *= sizeof(float);

   bytesUsed   += newLength - oldLength;
   bytesReused += oldLength;

   if (verboseAlloc)
      {
       FPR( stderr,
            "\n_MEM FLT_RE_AL Name: %s Free: %p",
            caller, ptr );
      }

   if (ptr)
      newPtr = (float *) realloc( (void *)ptr, (size_t)newLength );
   else
      newPtr = (float *) malloc( (size_t)newLength );

   if (verboseAlloc)
      {
       FPR( stderr, " Alloc: %p", newPtr );

#ifdef NMR64
       FPR( stderr, " Add: %ld Delete: %ld Now: %ld\n", newLength, oldLength, bytesUsed );
#else
       FPR( stderr, " Add: %d Delete: %d Now: %d\n", newLength, oldLength, bytesUsed );
#endif
      }

   ERRTEST( newPtr );

   return( newPtr );
}

/***/
/* Reallocate a type char array:
/***/

char *charReMalloc64( caller, ptr, newLength, oldLength )

   NMR_INT newLength, oldLength;
   char   *caller;
   char   *ptr;
{
   char *newPtr;

   if (!ptr) oldLength = 0;

   bytesUsed   += newLength - oldLength;
   bytesReused += oldLength;

   if (verboseAlloc)
      {
       FPR( stderr,
            "\n_MEM CHR_RE_AL Name: %s Free: %p",
            caller, ptr );
      }

   if (ptr)
      newPtr = (char *) realloc( (void *)ptr, (size_t)newLength );
   else
      newPtr = (char *) malloc( (size_t)newLength );

   if (verboseAlloc)
      {
       FPR( stderr, " Alloc: %p", newPtr );

#ifdef NMR64
       FPR( stderr, " Add: %ld Delete: %ld Now: %ld\n", newLength, oldLength, bytesUsed );
#else
       FPR( stderr, " Add: %d Delete: %d Now: %d\n", newLength, oldLength, bytesUsed );
#endif
      }

   ERRTEST( newPtr );

   return( newPtr );
}

/***/
/* Reallocate a type void array:
/***/

void *voidReMalloc64( caller, ptr, newLength, oldLength )

   NMR_INT newLength, oldLength;
   char    *caller;
   void    *ptr;
{
   void *newPtr;

   if (!ptr) oldLength = 0;

   bytesUsed   += newLength - oldLength;
   bytesReused += oldLength;

   if (verboseAlloc)
      {
       FPR( stderr,
            "\n_MEM CHR_RE_AL Name: %s Free: %p",
            caller, ptr );
      }

   if (ptr)
      newPtr = (void *) realloc( (void *)ptr, (size_t)newLength );
   else
      newPtr = (void *) malloc( (size_t)newLength );

   if (verboseAlloc)
      {
       FPR( stderr, " Alloc: %p", newPtr );

#ifdef NMR64
       FPR( stderr, " Add: %ld Delete: %ld Now: %ld\n", newLength, oldLength, bytesUsed );
#else
       FPR( stderr, " Add: %d Delete: %d Now: %d\n", newLength, oldLength, bytesUsed );
#endif
      }

   ERRTEST( newPtr );

   return( newPtr );
}

/***/
/* Deallocate argv-style Null-ptr terminated string list.
/***/

int freeStrList( s, n )

   char **s;
   int  n;
{
    int i;

    if (!s || n < 1) return( 0 );
    if (!*s) return( 0 );
 
    for( i = 0; i < n; i++ )
       {
        if (s[i])
           { 
            (void) strFree( s[i] );
            s[i] = (char *)NULL;
           }
       } 

    (void) deAlloc( "strDup", s, sizeof(char *)*(n+1) );

    return( 0 );
}

/***/
/* Allocate 2D pointer list given an allocated matrix.
/***/

float **Mat2Ptr64( matrix, xSize, ySize )

   NMR_INT xSize, ySize;
   float   *matrix;
{
    float **mList, **mPtr;

    mList = (float **) voidAlloc( "mat2ptr", sizeof(float *)*ySize );

    if (!mList)
       {
        return( mList );
       }

    mPtr = mList;

    while( ySize-- )
       {
        *mPtr++ = matrix;
        matrix += xSize;
       }
 
    return( mList ); 
}

/***/
/* Allocate 3D pointer list given an allocated matrix.
/***/

float ***Mat2Ptr3D64( matrix, xSize, ySize, zSize )

   NMR_INT xSize, ySize, zSize;
   float   *matrix;
{
    NMR_INT iy, iz;
    float   ***mList;

    mList = (float ***) voidAlloc( "mat2ptr", sizeof(float **)*zSize );

    if (!mList) return( mList );

    for( iz = 0; iz < zSize; iz++ )
       {
        mList[iz] = (float **) voidAlloc( "mat2ptr", sizeof(float *)*ySize );

        if (!mList[iz])
           {
            (void) deAlloc( "mat2ptr", mList, sizeof(float **)*zSize );
            return( (float ***)NULL );
           }

        for( iy = 0; iy < ySize; iy++ )
           {
            mList[iz][iy] = matrix;
            matrix        += xSize;
           }
       }

    return( mList );
}

/***/
/* Allocate 4D pointer list given an allocated matrix.
/***/

float ****Mat2Ptr4D64( matrix, xSize, ySize, zSize, aSize )

   NMR_INT  xSize, ySize, zSize, aSize;
   float   *matrix;
{
    NMR_INT iy, iz, ia;
    float   ****mList;

    mList = (float ****) voidAlloc( "mat2ptr", sizeof(float ***)*aSize );

    if (!mList) return( mList );

    for( ia = 0; ia < aSize; ia++ )
       {
        mList[ia] = (float ***) voidAlloc( "mat2ptr", sizeof(float **)*zSize );

        if (!mList[ia])
           {
            return( (float ****)NULL );
           }

        for( iz = 0; iz < zSize; iz++ )
           {
            mList[ia][iz] = (float **) voidAlloc( "mat2ptr", sizeof(float *)*ySize );

            if (!mList[ia][iz])
               {
                return( (float ****)NULL );
               }

            for( iy = 0; iy < ySize; iy++ )
               {
                mList[ia][iz][iy]  = matrix;
                matrix            += xSize;
               }
           }
       }

    return( mList );
}

/***/
/* Mat2PtrZ: generate pointer list for a sub-matrix.
/***/

float **Mat2PtrZ64( matrix, xSize, xOffset, yStart, yEnd )

   NMR_INT xSize, xOffset, yStart, yEnd; 
   float   *matrix;
{
    float   **mList, **mPtr;
    NMR_INT ySize;

    ySize = 1 + yEnd - yStart;

    mList = (float **) voidAlloc( "mat2ptr", sizeof(float *)*ySize );

    if (!mList)
       {
        return( mList );
       }

    matrix += (yStart - 1)*xSize + xOffset - 1;
    mPtr    = mList;

    while( ySize-- )
       {
        *mPtr++ = matrix;
        matrix += xSize;
       }

    return( mList );
}

/***/
/* Allocate and return a 2D pointer array describing
/* an already-allocated 2D double precision matrix.
/***/

double **Mat2PtrD64( matrix, xSize, ySize )

   NMR_INT xSize, ySize;
   double  *matrix;
{
    double **mList, **mPtr;

    mList = (double **) voidAlloc( "mat2ptr", sizeof(double *)*ySize );

    if (!mList)
       {
        return( mList );
       }

    mPtr = mList;

    while( ySize-- )
       {
        *mPtr++ = matrix;
        matrix += xSize;
       }

    return( mList );
}

/***/
/* Deallocate a 2D matrix and associated pointers.
/***/

int MatFree64( matrixPtr, xSize, ySize )

   NMR_INT xSize, ySize;
   float   **matrixPtr;
{
    if (matrixPtr)
       {
        (void) deAlloc( "matFree", *matrixPtr, sizeof(float)*xSize*ySize );
        (void) deAlloc( "matFree", matrixPtr, sizeof(float *)*ySize );
       }

    return( 0 );
}

/***/
/* Deallocate pointers associated with a 3D matrix.
/***/

int MatPtrFree3D64( matrixPtr, xSize, ySize, zSize )

   NMR_INT xSize, ySize, zSize;
   float   ***matrixPtr;
{
    NMR_INT iz;

    if (!matrixPtr) return( 0 );

    for( iz = 0; iz < zSize; iz++ )
       {
        (void) deAlloc( "mat2ptr", matrixPtr[iz], sizeof(float *)*ySize );
       }

    (void) deAlloc( "mat2ptr", matrixPtr, sizeof(float **)*zSize );
 
    return( 0 );
}

/***/
/* Deallocates pointers associated with a 4D matrix.
/***/

int MatPtrFree4D64( matrixPtr, xSize, ySize, zSize, aSize )

   NMR_INT xSize, ySize, zSize, aSize;
   float   ****matrixPtr;
{
    NMR_INT iz, ia;

    if (!matrixPtr) return( 0 );

    for( ia = 0; ia < aSize; ia++ )
       {
        if (!matrixPtr[ia]) break;

        for( iz = 0; iz < zSize; iz++ )
           {
            if (!matrixPtr[ia][iz]) break;

            (void) deAlloc( "mat2ptr", matrixPtr[ia][iz], sizeof(float *)*ySize );
           }

        (void) deAlloc( "mat2ptr", matrixPtr[ia], sizeof(float **)*zSize );
       }

    (void) deAlloc( "mat2ptr", matrixPtr, sizeof(float ***)*aSize );

    return( 0 );
}

/***/
/* Deallocate a 2D double prec matrix and associated pointers.
/***/

int MatFreeD64( matrixPtr, xSize, ySize )

   NMR_INT xSize, ySize;
   double  **matrixPtr;
{
    if (matrixPtr)
       {
        (void) deAlloc( "matFree", *matrixPtr, sizeof(double)*xSize*ySize );
        (void) deAlloc( "matFree", matrixPtr,  sizeof(double *)*ySize );
       }

    return( 0 );
}

/***/
/* Deallocate memory without updating debug tallies of memory used.
/***/

void unAllocS64( caller, ptr, length )

   NMR_INT length;
   char    *caller;
   void    *ptr;
{
   if (!ptr || length < 1) return;
   (void) unAlloc64( caller, ptr, length );
}

/***/
/* Deallocate memory.
/***/

void unAlloc64( caller, ptr, length )

   NMR_INT length;
   char    *caller;
   void    *ptr;
{
   if (ptr) bytesUsed -= length;

   if (verboseAlloc)
      {
#ifdef NMR64
       FPR( stderr,
            "\n_MEM %s Name: %s Free: %p Delete: %ld Now: %ld Reused: %ld\n",
            ptr ? "DEALLOC" : "NOOP_AL", caller, ptr,
            length, bytesUsed, bytesReused );
#else
       FPR( stderr,
            "\n_MEM %s Name: %s Free: %p Delete: %d Now: %d Reused: %d\n",
            ptr ? "DEALLOC" : "NOOP_AL", caller, ptr,
            length, bytesUsed, bytesReused );
#endif
      }

   if (ptr) (void) free( ptr );
}

/***/
/* Bottom.
/***/
