/***/
/* atof.h: substitute for atof, which fails on some systems
/*         given integer input text.
/***/

#include <stdlib.h>

#include "prec.h"

#ifndef __atof_h
#define __atof_h

#if defined (CONVEX) || defined (LINUX) || defined(WINXP)
#define ATOF Atof2
double Atof2();
#else
#define ATOF atof
#endif

#define ATOI atoi

#ifdef NMR64
#define ATON atol
#else
#define ATON atoi
#endif
 
#endif
