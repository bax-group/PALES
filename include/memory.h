/***/
/* Definition of dynamic memory wrappers:
/***/

#include "prec.h"

#ifndef __nmr_memory_h

#define __nmr_memory_h

#ifdef LINUX
   void free();
#endif

#ifdef WINXP 
   void free();
#endif

   int    setAlloc();
   int    getMemUsed();
   int    showAlloc();
   int    ShowCurrentAlloc();

   void     *voidMalloc64();
   char     *charMalloc64();
   float    *fltMalloc64();
   int      *intMalloc64();
   double   *dFltMalloc64();
   long int *dIntMalloc64();
   NMR_INT  *n64Malloc();

   void  *voidReMalloc64();
   char  *charReMalloc64();
   float *fltReMalloc64();
   int   *intReMalloc64();

   char  *strDup();
   char  *strDup2();
   int   strFree();
   int   freeStrList();

   float  **MatMalloc64();
   double **MatMallocD64();

   double **Mat2PtrD64();
   float  **Mat2PtrZ64();

   float  **Mat2Ptr64();
   float  ***Mat2Ptr3D64();
   float  ****Mat2Ptr4D64();

   int    MatPtrFree3D64();
   int    MatPtrFree4D64();

   int    MatFree64();
   int    MatFreeD64();

   void   unAlloc64();
   void   unAllocS64();

#ifdef VOID
#undef VOID
#endif

#ifdef SOLARIS
#define VOID   void
#define IOPTR  void 
#define IOSIZE size_t
#define SIZE_T size_t 
#define SIZE_DEFS
#endif

#if defined (WINNT) || defined (WIN95)
#define VOID   void
#define IOPTR  void 
#define IOSIZE size_t
#define SIZE_T size_t 
#define SIZE_DEFS
#endif

#ifdef SUN
#ifndef SIZE_DEFS
#define VOID   void
#define IOPTR  char
#define IOSIZE int
#define SIZE_T unsigned
#define SIZE_DEFS
#endif
#endif

#ifdef SGI
#define VOID   void
#define IOPTR  void 
#define IOSIZE unsigned
#define SIZE_T size_t 
#define SIZE_DEFS
#endif

#ifdef IBM 
#define VOID   void
#define IOPTR  void
#define IOSIZE size_t
#define SIZE_T size_t
#define SIZE_DEFS
#endif

#ifdef HP 
#define VOID   void
#define IOPTR  void
#define IOSIZE size_t
#define SIZE_T size_t
#define SIZE_DEFS
#endif

#ifdef ALPHA
#define VOID   void
#define IOPTR  char
#define IOSIZE size_t
#define SIZE_T size_t
#define SIZE_DEFS
#endif

#ifndef SIZE_DEFS
#define VOID   void
#define IOPTR  void 
#define IOSIZE size_t
#define SIZE_T size_t 
#define SIZE_DEFS
#endif

#define voidAlloc( WHO, LEN )  voidMalloc64( WHO, (NMR_INT)((NMR_INT)LEN) )
#define charAlloc( WHO, LEN )  charMalloc64( WHO, (NMR_INT)((NMR_INT)LEN) )
#define fltAlloc( WHO, LEN )   fltMalloc64(  WHO, (NMR_INT)((NMR_INT)LEN) )
#define intAlloc( WHO, LEN )   intMalloc64(  WHO, (NMR_INT)((NMR_INT)LEN) )
#define dFltAlloc( WHO, LEN )  dFltMalloc64( WHO, (NMR_INT)((NMR_INT)LEN) )
#define dIntAlloc( WHO, LEN )  dIntMalloc64( WHO, (NMR_INT)((NMR_INT)LEN) )
#define n64Alloc( WHO, LEN )   n64Malloc(  WHO, (NMR_INT)((NMR_INT)LEN) )

#define voidReAlloc(WHO,PTR,NEWLEN,OLDLEN) voidReMalloc64( WHO, (void *)PTR,  (NMR_INT)((NMR_INT)NEWLEN), (NMR_INT)((NMR_INT)OLDLEN) )
#define charReAlloc(WHO,PTR,NEWLEN,OLDLEN) charReMalloc64( WHO, (char *)PTR,  (NMR_INT)((NMR_INT)NEWLEN), (NMR_INT)((NMR_INT)OLDLEN) )
#define fltReAlloc(WHO,PTR,NEWLEN,OLDLEN)  fltReMalloc64(  WHO, (float *)PTR, (NMR_INT)((NMR_INT)NEWLEN), (NMR_INT)((NMR_INT)OLDLEN) )
#define intReAlloc(WHO,PTR,NEWLEN,OLDLEN)  intReMalloc64(  WHO, (int *)PTR,   (NMR_INT)((NMR_INT)NEWLEN), (NMR_INT)((NMR_INT)OLDLEN) )

#define deAlloc(WHO,PTR,LEN)      unAlloc64(  WHO, (void *)PTR, (NMR_INT)((NMR_INT)LEN) )
#define deAllocS(WHO,PTR,LEN)     unAllocS64( WHO, (void *)PTR, (NMR_INT)((NMR_INT)LEN) )

#define matAlloc( WHO, NX, NY )        MatMalloc64(  WHO, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY) )
#define matAllocD( WHO, NX, NY )       MatMallocD64( WHO, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY) )

#define mat2Ptr( M, NX, NY )           Mat2Ptr64(   M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY) )
#define mat2Ptr3D( M, NX, NY, NZ )     Mat2Ptr3D64( M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY), (NMR_INT)((NMR_INT)NZ) )
#define mat2Ptr4D( M, NX, NY, NZ, NA ) Mat2Ptr4D64( M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY), (NMR_INT)((NMR_INT)NZ), (NMR_INT)((NMR_INT)NA) )

#define mat2PtrD( M, NX, NY )           Mat2PtrD64( M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY) )
#define mat2PtrZ( M, NX, XOFF, Y1, YN ) Mat2PtrZ64( M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)XOFF), (NMR_INT)((NMR_INT)Y1), (NMR_INT)((NMR_INT)YN) )

#define matFree( M, NX, NY )              MatFree64( M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY) )

#define matPtrFree3D( M, NX, NY, NZ )     MatPtrFree3D64( M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY), (NMR_INT)((NMR_INT)NZ) )
#define matPtrFree4D( M, NX, NY, NZ, NA ) MatPtrFree4D64( M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY), (NMR_INT)((NMR_INT)NZ), (NMR_INT)((NMR_INT)NA) )

#define matFreeD( M, NX, NY )          MatFreeD64( M, (NMR_INT)((NMR_INT)NX), (NMR_INT)((NMR_INT)NY) )

#define showCurrentAlloc( TYPE, WHO, PTR, N ) ShowCurrentAlloc64( TYPE, WHO, (void *)PTR, (NMR_INT)((NMR_INT)N) )

#define MTEST( TXT, EXP )                                                 \
   if (!(EXP))                                                            \
      {                                                                   \
       (void) fprintf( stderr, "%s Error allocating memory.\n", TXT );    \
       error = 1;                                                         \
       goto shutdown;                                                     \
      } error = 0

#endif /* __nmr_memory_h */
