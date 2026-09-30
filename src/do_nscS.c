#include "nr.h"
#include "nsc.h"
#include <stdlib.h>
#include <stddef.h>

#include "cmndargs.h"
#include "memory.h"
#include "prec.h"

#include "hsymbol.h"
#include "pdb.h"
#include "dc.h"

#define FPR (void)fprintf

#define SAT (0.5)		/* 0.57 */
#define VDW (1.4)
#define radiusUni (1.7)


int surfSel( int atomCount, double *xI, double *yI,
		double *zI, int *atomSelCount );



int surfSel( int atomCount, double *xI, double *yI,
		double *zI, int *atomSelCount )
{
  int      error;
  int      i, n, k, density = 700, n_dots;
  double   area, volume, *atomwise_area, *surface_dots;
  double   *coordinates, *radii;
  double   *xJ, *yJ, *zJ;


/***/
/* Initialize:
/***/
   error   = 0;

   i       = 0;
   k       = 0;
   n       = 0;
   density = 700;
   n_dots  = 0;
   area    = 0.0;
   volume  = 0.0;

   xJ            = (double *)NULL;
   yJ            = (double *)NULL;
   zJ            = (double *)NULL;
   radii         = (double *)NULL;
   coordinates   = (double *)NULL;
   atomwise_area = (double *)NULL;
   surface_dots  = (double *)NULL;


/***/
/* Allocate memory for coordinates and radii
/***/

    radii          = dFltAlloc( "surfSel", atomCount );
    coordinates    = dFltAlloc( "surfSel", atomCount*3 );
    xJ             = dFltAlloc( "surfSel", atomCount+1 );
    yJ             = dFltAlloc( "surfSel", atomCount+1 );
    zJ             = dFltAlloc( "surfSel", atomCount+1 );
 
    if ( !radii || !coordinates ||!xJ ||!yJ ||!zJ )
       {
        FPR( stderr, "surfSel: Error allocating coord buffers.\n" );
        error = 1;
        goto shutdown;
      }


/***/
/* Store coordinates and radii temporarily
/* Record index of stored values
/***/
    for( i = 0; i < atomCount; i++ )
       {
	   xJ[i+1] = xI[i+1];
	   yJ[i+1] = yI[i+1];
	   zJ[i+1] = zI[i+1];
	   coordinates[i*3]   = xI[i+1];
	   coordinates[i*3+1] = yI[i+1];
	   coordinates[i*3+2] = zI[i+1];
           radii[i] = radiusUni + VDW;
       }



/***/
/* Calculate dot surface
/***/
    if (NSC(coordinates, radii, atomCount, density,
          FLAG_ATOM_AREA | FLAG_VOLUME,
          &area, &atomwise_area, &volume, &surface_dots, &n_dots))
       {
          FPR( stderr, "NSC: error.\n" );
          error = 1;
          goto shutdown;
       }



/***/
/* Select surface accessible atoms
/***/

  k = 0;
  for( i = 0; i < atomCount; i++ )
    {
	if(atomwise_area[i]/(4*PI*(radii[i]-VDW)*(radii[i]-VDW))>SAT)
	  {
	     k++;
	     xI[k] = xJ[i+1];
	     yI[k] = yJ[i+1];
	     zI[k] = zJ[i+1];
	  }
    }
  *atomSelCount = k;



/***/
/* Special cases:
/***/

shutdown:

   (void) deAlloc( "surfSel", xJ,          sizeof(double)*atomCount+1 );
   (void) deAlloc( "surfSel", yJ,          sizeof(double)*atomCount+1 );
   (void) deAlloc( "surfSel", zJ,          sizeof(double)*atomCount+1 );
   (void) deAlloc( "surfSel", radii,       sizeof(double)*atomCount );
   (void) deAlloc( "surfSel", coordinates, sizeof(double)*atomCount*3 );

   return( error );
}
/* ========================================================================= */

#undef SAT
#undef PI
#undef VDW
#undef radiusUni
