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


int callNSC( struct PDBInfo *pdbInfo );
float getRadius( struct AtomInfo *atomInfo );



int callNSC( struct PDBInfo *pdbInfo )
{
  int      atomCount, selectCount, error;
  int      i, n, k, density = 700, n_dots;
  double   area, volume, *atomwise_area, *surface_dots;
  int      *indAtomList;
  double   *coordinates, *radii;
  struct   AtomInfo *atomList;


/***/
/* Initialize:
/***/
   error   = 0;

   i       = 0;
   n       = 0;
   k       = 0;
   density = 700;
   n_dots  = 0;
   area    = 0.0;
   volume  = 0.0;

   atomCount = 0;
   selectCount = 0;
   atomList  = (struct AtomInfo *)NULL;

   indAtomList   = (int    *)NULL;
   radii         = (double *)NULL;
   coordinates   = (double *)NULL;
   atomwise_area = (double *)NULL;
   surface_dots  = (double *)NULL;


/***/
/* Allocate memory for coordinates and radii
/***/

    atomCount = pdbInfo->atomCount;
    atomList  = pdbInfo->atomList;

    indAtomList    = intAlloc( "callNSC", atomCount );
    radii          = dFltAlloc( "callNSC", atomCount );
    coordinates    = dFltAlloc( "callNSC", atomCount*3 );
 
    if (!radii || !coordinates|| !indAtomList)
       {
        FPR( stderr, "callNSC: Error allocating coord buffers.\n" );
        error = 1;
        goto shutdown;
      }


/***/
/* Store coordinates and radii temporarily
/* Record index of stored values
/***/
    k = 0;
    for( i = 0; i < atomCount; i++ )
       {
	if(atomList[i].selected)
	   {
	   coordinates[k*3]   = atomList[i].x;
	   coordinates[k*3+1] = atomList[i].y;
	   coordinates[k*3+2] = atomList[i].z;
           radii[k] = getRadius( atomList + i ) + VDW;
 	   indAtomList[k] = i;
	   k++;
	   }
       }



/***/
/* Calculate dot surface
/***/
    if (NSC(coordinates, radii, k, density,
          FLAG_ATOM_AREA | FLAG_VOLUME,
          &area, &atomwise_area, &volume, &surface_dots, &n_dots)) {
        FPR( stderr, "NSC: error.\n" );
        error = 1;
        goto shutdown;
  }



/***/
/* Select surface accessible atoms
/***/

  selectCount = pdbSetSelection( pdbInfo, 0 );

  for( i = 0; i < k; i++ )
    {
	n = indAtomList[i];

	if(atomwise_area[i]/(4*PI*(radii[i]-VDW)*(radii[i]-VDW))>SAT) {
		selectCount = pdbSelect1( pdbInfo,
			atomList[n].resID, 
			atomList[n].resID, 
			atomList[n].atomName, 1, 1 );

	}
    }


/***/
/* Deallocate memory
/***/
   (void) deAlloc( "callNSC", indAtomList, sizeof(int)*atomCount );
   (void) deAlloc( "callNSC", radii,       sizeof(double)*atomCount );
   (void) deAlloc( "callNSC", coordinates, sizeof(double)*atomCount*3 );


/***/
/* Special cases:
/***/

shutdown:

   return( error );
}
/* ========================================================================= */


/* ========================================================================= */
float getRadius( struct AtomInfo *atomInfo )
{
    float r;

    switch( atomInfo->atomName[0] )
       {
        case 'C':
           r = 1.9;
           break;
        case 'F':
           r = 1.9;
           break;
        case 'H':
           r = 1.4;
           break;
        case 'N':
           r = 1.8;
           break;
        case 'O':
           r = 1.4;
           break;
        case 'P':
           r = 1.8;
           break;
        case 'S':
           r = 1.8;
           break;
        default:
           r = 1.4;
           break;
       }

    return( r );
}
/* ========================================================================= */

#undef SAT
#undef PI
#undef VDW
