#ifndef _nmr_prec_h
#define _nmr_prec_h

#define INT4  int
#define REAL4 float
#define REAL8 double
#define CHAR  char

#ifdef NMR64
#define INT8         long
#define UINT8        unsigned long
#define NMR_INT      long
#define NMR_INT_FMT  %ld
#define A_TO_NMR_INT atol
#else
#define INT8         long long
#define UINT8        unsigned long long
#define NMR_INT      int
#define NMR_INT_FMT  %d
#define A_TO_NMR_INT atoi
#endif

#define NAMELEN          1023 
#define BYTESPERWORD        4
#define BITSPERWORD        32
#define LARGENUM     10.0E+16

#ifdef NMR64
#define NMR_INT_MAX SSIZE_MAX
#else
#define NMR_INT_MAX INT_MAX
#endif

#ifdef PI
#undef PI
#endif

#define PI 3.14159265

#ifdef WIN95
#define strcasecmp(S1,S2) strcmpi(S1,S2)
#endif

#undef SYS_BYTESWAP

#ifdef PC_SOLARIS
#define SYS_BYTESWAP
#endif

#ifdef WINNT
#define SYS_BYTESWAP
#endif

#ifdef WIN95
#define SYS_BYTESWAP
#endif

#ifdef LINUX 
#define SYS_BYTESWAP
#endif

#ifdef WINXP 
#define SYS_BYTESWAP
#endif

#ifdef ALPHA 
#define SYS_BYTESWAP
#endif

#ifdef MAC_OSX
#ifdef __LITTLE_ENDIAN__
#define SYS_BYTESWAP
#endif
#endif

#endif
