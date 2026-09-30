/*********************************************************/
/*********************************************************/
/* pales.c simulates the alignment of macromolecules     */
/* based on their three-dimensional structure            */
/* (C) Markus Zweckstetter                 Version 2.1   */
/*********************************************************/
/*********************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "pales.h"
#include "unsp.h"

#include "cmndargs.h"
#include "memory.h"
#include "prec.h"
#include "atof.h"

#include "stralloc.h"
#include "strmatch.h"
#include "hsymbol.h"
#include "pdb.h"
#include "pdbtab.h"
#include "dc.h"

#include "dynamo.h"

#define FPR             (void)fprintf
#define SPR             (void)sprintf
#define ISWAP( IA, IB ) itemp = IA; IA = IB; IB = itemp

static int dbgFlag    = 0;

#define RAD(x) (x*PI/180.)
#define DEG(x) (x*180./PI)
#define S(x) sin(x)
#define C(x) cos(x)


#define epsilon		7.82853e+01  /* unitless */
#define k_B		8.31451e-03  /* [kJ/mol.K] */
#define e_2		1.38935e+02  /* [nm.kJ/mol] */

#define gH        2.675198e+8
#define gN       -2.7116e+7
#define hbCons    1.05457266e-4		/* as given in Bronstein */
#define u0	  PI*4.0e-7
#define cx	  ((u0*hbCons)/(4.0*PI*PI))
#define diNH (-1*gH*gN*cx/(RAD_N_HN*RAD_N_HN*RAD_N_HN))

float Pmax=10.0;


/* ========================================================================= */
/* ========================================================================= */
/* ========================================================================= */


/* ======================== rotate vector ================================= */
void rotate_vector(double theta, double phi, double psi, double *input_vector,
        double *rotated_vector);


/* ================== rotate coordinate system =========================== */
void rotate_coorsystem(double theta, double phi, double psi, double *x_axis_tensor,
        double *y_axis_tensor, double *z_axis_tensor);


/* ======= calculation of Szz according to Prestegard ===================== */
double Szz(double tensor_zaxis[4]);


/* ========== calculation of Sxx according to Prestegard =============== */
double Sxx(double tensor_xaxis[4]);


/* ========== calculation of Syy according to Prestegard =============== */
double Syy(double tensor_yaxis[4]);


/* ============= calculation of Sxy according to Prestegard ============== */
double Sxy(double tensor_xaxis[4], double tensor_yaxis[4]);


/* ============ calculation of Sxz according to Prestegard =============== */
double Sxz(double tensor_xaxis[4], double tensor_zaxis[4]);


/* ============ calculation of Syz according to Prestegard =============== */
double Syz(double tensor_yaxis[4], double tensor_zaxis[4]);


/* ============== simulate sterically induced alignment =================== */
void ssia( struct PalesInfo *palesInfo, int no_atoms, 
           double *x, double *y, double *z,
           double max_protein_length, double min_protein_length,
           int n_dots, double *dots, double *saupe );

int   refine_bicelle( struct PalesInfo *palesInfo );

int   make_eigVecList(double eigVec[3][3], double eigVecXlist[NO_LIST][3],
        double eigVecYlist[NO_LIST][3], double eigVecZlist[NO_LIST][3]);

int   cross_productM2(int n, int m, double eigVec[3][3], double product[3]);

int   add_eigVec_List(double eigVecXlist[NO_LIST][3], double eigVecYlist[NO_LIST][3],
        double eigVecZlist[NO_LIST][3], int n, double eigVec[3][3]);

int   mirror_solutions(int n, double eigVec[3][3], double eigVecXlist[NO_LIST][3],
        double eigVecYlist[NO_LIST][3], double eigVecZlist[NO_LIST][3]);

void  legalNotice();

/* ========================================================================= */
/* ========================================================================= */

int   getDCTensorParams( int argc, char **argv, struct DTensor *tInfo );

int   printInfo( );

int   coorSelect( struct PDBInfo *pdbInfo, struct PalesInfo *palesInfo,
                int *atomSelCount, 
                double *xI, double *yI, double *zI );

int   chCoorSelect(struct PDBInfo *pdbInfo, struct PalesInfo *palesInfo, 
               int *atomSelCount, 
	       double *xI, double *yI, double *zI ,
               int chCount, int *chargeAtomID);

int   simulateAlignment( struct PalesInfo *palesInfo, int atomSelCount, 
        double *xI, double *yI, double *zI , 
	int chCount, int *chargeAtomID, float *charges,
	int dotCount, double *ddots,
	int digPotential, double *potential, double offset, double stepsize,
	double saupe[5] );

int   readPDBfree(char *inName, int *atomCount, double *xI, 
	double *yI, double *zI);

int   rotateDomain( double angle[3], struct PDBInfo *pdbInfo,
	int resN, int resC );

int    positionDomain(struct PDBInfo *pdbInfo, struct PalesInfo *palesInfo,
        int resN, int resC,
        int origin1, int origin2 );

int    writeEulerInfo( double irrSaupe[5], double gMag,
	double eigVal[3], double eigVec[3][3],
	double alphaList[NO_LIST], double betaList[NO_LIST],
	double gammaList[NO_LIST], char *result, int verbose );

int    calcDcos(double eigVec [3][3], double eigVec_svd[3][3], 
	float dcos[3]);

int    simDC( struct PDBInfo *pdbInfo, struct DTable *dTab,
           struct PalesInfo *palesInfo, char *remark,
           char *outDName, char *pdbRotName, int verbose );

int    mapSimDC( struct PDBInfo *pdbInfo, struct DTable *dTab, 
           struct PalesInfo *palesInfo, 
           char *outDName, char *pdbRotName, int verbose,
           char *outMapName, int mapCount, char *remark);

int    selectSurfDom( struct PDBInfo *pdbInfo, int resN, int resC,
	 int *indICount, int *indI, int *indJCount, int *indJ );

float  calcNorm( float x,  float y,  float z );

int    ssiaFree( struct PalesInfo *palesInfo, char *pdbFreeName,
                char *outAName, struct DTable *dTab );

int    shiftMol2Domain( struct PDBInfo *pdbInfo, int resN, int resC );

double gauss(void);
double genrand(void);

int    compareSaupe( struct PalesInfo *palesInfo,
                double userSaupe1[5], double userSaupe2[5],
		char *outName, struct DTable *dTab, char *remark);

int    addSaupe( int opCode, float c1, float c2,
                double userSaupe1[5], double userSaupe2[5],
                double saupe[5] );

int    selTypDC( struct DInfo *dList, int dCount, char *atomNameI, 
                  char *atomNameJ );

int    analyzeDC( char *inD1Name, char *inD2Name,
               struct PalesInfo *palesInfo, struct DcInfo *dcInfo,
	       char *atomNameS1, char *atomNameS2, char *outName);

int    adjustDC( int opCode, float c1, float c2,
              struct DTable *dTab1, struct DTable *dTab2 );

float  addDC( int opCode, float c1, float c2,
                float d1, float d2 );

int    palesInfoNull( struct PalesInfo *palesInfo );

int    dcInfoNull( struct DcInfo *dcInfo );

int    gridInfoNull( struct GridInfo *gridInfo );

int    elInfoNull( struct ElInfo *elInfo );
 
int    writePales( struct PalesInfo *palesInfo, char *result );

int    getStatDCmod(struct DTable *dTab, 
          float *rSumD2, float *pSumD, float *qD, float *coorR,
          float *offset,  float *errOffset, float *slope, float *errSlope,
          float *aveSlope, float *errAveSlope,
          int *nD );

int    getEulerDC (float eigenVec[3][3], struct DTensor *dTensor);

/* ========================================================================= */
/* ========================================================================= */


int main(int argc, char *argv[])
{
  int i;

  /* UNSP */
  int    density = 100, n_dots=0, n_ddots=0;
  double *dots, *ddots;

  /* electrostatics */
  int    digPotential;            /* digitization of potential */
  double offset,stepsize;         /* starting position and grid size of potential */
  int    *chargeAtomID;
  float  *charges;
  double *potential;
 
  /* analyzation of Saupe matrix */
  static double saupe[5], irrSaupe[5], saupeCp1[5], saupeCp2[5];
  static double eigVal[3],eigVec[3][3];
  static double gMag,da,dr,R;
  static double alphaList[NO_LIST], betaList[NO_LIST], gammaList[NO_LIST];

  double  rotAngle[3];


  struct  PalesInfo  *palesInfo;
  struct  DcInfo     *dcInfo;


/*** FRANK ***/
    struct AtomInfo *atomList;
    struct PDBInfo  *pdbInfo;
    struct DTable   *dTab;
    struct CSTable  *csTab;


    float  rSum, pSum;
    double *xI, *yI, *zI;
    char   inpName[NAMELEN+1];
    char   potName[NAMELEN+1], pkaName[NAMELEN+1];
    char   inEName[NAMELEN+1], pdbFreeName[NAMELEN+1]; 
    char   outEAxName[NAMELEN+1], flexName[NAMELEN+1];
    char   pdbRotName[NAMELEN+1], outDcEnsName[NAMELEN+1];
    char   aveName[NAMELEN+1], outEName[NAMELEN+1];
    char   outAName[NAMELEN+1];

    int    dcCount, atomCount, atomSelCount, chCount;

    int    validCount, eulerCount, aOff, n;
    int    mapCount;

    char   result[8000], ctemp[8000], *remark;

    int    dCount, csCount, validCountD, validCountCS;
    int    fitFlag, wCount, itemp, error;

    char   inDName[NAMELEN+1], inCSName[NAMELEN+1], outDName[NAMELEN+1],
           outCSName[NAMELEN+1], outSName[NAMELEN+1], pdbList[NAMELEN+1],
           pdbPath[NAMELEN+1], pdbName[NAMELEN+1], refName[NAMELEN+1],
           surfPath[NAMELEN+1], atomNames[NAMELEN+1], resNames[NAMELEN+1],
           rcName[NAMELEN+1], homoName[NAMELEN+1], outMapName[NAMELEN+1],
           atomNameS1[NAMELEN+1], atomNameS2[NAMELEN+1],
           inD1Name[NAMELEN+1], inD2Name[NAMELEN+1];

    float  *wList, csThresh, validThresh;

    struct DInfo   *dObsList, *dCalcList;
    struct  PDBSelectInfo sInfo;

int selectCount;

/***/
/* Initialize:
/***/

    pdbInfo   = (struct PDBInfo *)NULL;
    atomList  = (struct AtomInfo *)NULL;
    dTab      = (struct DTable *)NULL;
    csTab     = (struct CSTable *)NULL;
    palesInfo = (struct PalesInfo *)NULL;
    dcInfo    = (struct DcInfo *)NULL;

    error     = 0;
    mapCount  = 0;
    selectCount = 0;
 
    (void) strcpy( pdbList,    "pdb.list" );
    (void) strcpy( rcName,     "rc.tab" );
    (void) strcpy( homoName,   "homo.tab" );
    (void) strcpy( atomNames,  "All" );
    (void) strcpy( resNames,   "None" );
    (void) strcpy( atomNameS1, "All" );
    (void) strcpy( atomNameS2, "All" );

    *surfPath     = '\0';
    *pdbPath      = '\0';
    *pdbName      = '\0';
    *refName      = '\0';
    *inDName      = '\0';
    *inCSName     = '\0';
    *outDName     = '\0';
    *outCSName    = '\0';
    *outSName     = '\0';

    *outMapName   = '\0';
    *inpName      = '\0';
    *potName      = '\0';
    *pkaName      = '\0';
    *inEName      = '\0';
    *aveName      = '\0';
    *flexName     = '\0';
    *pdbFreeName  = '\0';
    *pdbRotName   = '\0';
    *outEName     = '\0';
    *outEAxName   = '\0';
    *outDcEnsName = '\0';

    *outAName   = '\0';
    *inD1Name   = '\0';
    *inD2Name   = '\0';


    csThresh    = 2.00;
    validThresh = 0.95;

    fitFlag     = 0;

    wCount    = 0;
    wList     = (float *)NULL;

    (void) pdbSelectNull( &sInfo );

    (void) surfNull();
    (void) dipolarNull();


    *result   = '\0';
    *ctemp    = '\0';
    remark   = (char *)NULL;

    dcCount       = 0;
    atomCount     = 0;
    atomSelCount  = 0;
    chCount       = 0;

  digPotential = 0;
  offset       = 0.0;
  stepsize     = 0.0;



    xI           = (double *)NULL;
    yI           = (double *)NULL;
    zI           = (double *)NULL;
    chargeAtomID = (int *)NULL;
    charges      = (float *)NULL;
    potential    = (double *)NULL;
    dots         = (double *)NULL;
    ddots        = (double *)NULL;


    (void) dynInit( DYN_IGNORE_KEY, DYN_IGNORE_RAMA );


/***/
/* Program start:
/***/

  legalNotice();



/***/
/* Allocate PDB and DC structures.
/* Null-initialize them.
/***/

    dTab      = (struct DTable *)voidAlloc( "matrix",  sizeof(struct DTable) );
    csTab     = (struct CSTable *)voidAlloc( "matrix",  sizeof(struct CSTable) );
    pdbInfo   = (struct PDBInfo *)voidAlloc( "matrix", sizeof(struct PDBInfo) );
    palesInfo = (struct PalesInfo *)voidAlloc( "matrix",  sizeof(struct PalesInfo) );
    dcInfo    = (struct DcInfo *)voidAlloc( "matrix",  sizeof(struct DcInfo) );

    if (!dTab || !pdbInfo || !csTab || !palesInfo || !dcInfo)
       {
        FPR( stderr, "Error allocating DC/PDB/PALES structures.\n" );
        error = 1;
        goto shutdown;
       }

    (void) dTabNull( dTab );
    (void) csTabNull( csTab );
    (void) pdbNull( pdbInfo );
    (void) palesInfoNull( palesInfo );
    (void) dcInfoNull( dcInfo );


/***/
/* Extract command line:
/***/

    (void) initArgs( argc, argv, CHECK_ENV );


    error = getparams( argc, argv,
                       inpName, potName, pkaName, inEName, flexName,
                       outDcEnsName, pdbFreeName, pdbRotName, outEName,
                       outEAxName, aveName,
                       inDName, inCSName, outDName, outCSName, outSName,
                       surfPath, atomNames, resNames, pdbList, pdbPath,
                       pdbName, refName, &dTab->dTensor, &fitFlag,
                       &wList, &wCount, &csThresh, &validThresh,
                       &mapCount, outMapName,
		       saupeCp1, saupeCp2, outAName, &remark,
                       atomNameS1, atomNameS2, inD1Name, inD2Name,
                       palesInfo, &sInfo, dcInfo);       



    if (error)
       {
        FPR( stderr, "DC Error extracting command-line arguments.\n" );
        goto shutdown;
       }


    (void) checkArgs( argc, argv );

    if (dbgFlag > 1) (void) setAlloc( 1 );



/***/
/***/
/***/
/* Selection of mode:
/***/
/***/
/***/

    if ( (palesInfo->palesMode == PALES_STERIC_FREE) ||
         (palesInfo->palesMode == PALES_SAUPE_CP) ||
         (palesInfo->palesMode == PALES_DC_CP) )
      {


         switch( palesInfo->palesMode )
            {

/***/
/***/
/* Simulation of alignment if PDB file in free format
/***/
/***/

             case PALES_STERIC_FREE:

               if( error = ssiaFree( palesInfo, pdbFreeName, outAName, dTab ) )
                  {
                     FPR( stderr, "matrix: error in ssiaFree.\n" );
                     goto shutdown;
                  }

                break;
 


/***/
/***/
/* Compare 2 Saupe matrices :
/***/
/***/

             case PALES_SAUPE_CP:

               if( error = compareSaupe(palesInfo, saupeCp1, saupeCp2,
		   outAName, dTab, remark) )
                 {
                    FPR( stderr, "Matrix: error in comparison of order matrices.\n");
                    goto shutdown;
                 }

                break;



/***/
/***/
/* Compare 2 sets of dipolar couplings :
/***/
/***/

             case PALES_DC_CP:

               if( error = analyzeDC( inD1Name, inD2Name, palesInfo, dcInfo,
		   atomNameS1, atomNameS2, outDName) )
                 {
                    FPR( stderr, "Matrix: error in comparison of order matrices.\n");
                    goto shutdown;
                 }

               break;

             default:
               break;


            } /* end of switch PALES mode */


     }	/* end of 'if' unformatted functions */



/***/
/***/
/* Functions using formatted PDB file :
/***/
/***/


   else
     {

/***/
/* Read in PDB file:
/* Read in DC table:
/***/


        if( *pdbName )
          {
            if (error = readPDB( pdbName, pdbInfo, PDB_SELECT_NULL ))
               {
                 FPR( stderr, "Error reading PDB input %s\n", pdbName );
                 goto shutdown;
               }

            FPR( stderr,
                  "%s, %d Atoms, Residues %d to %d.\n",
                  pdbName, pdbInfo->atomCount, pdbInfo->minResID, pdbInfo->maxResID );
          }



        if (*inDName)
           {
            if (error = readDC( inDName, dTab ))
               {
                  FPR( stderr, "Error reading DC input %s\n", inDName );
                  goto shutdown;
               }

            FPR( stderr,
                  "%s, %d Couplings, Residues %d to %d.\n",
                  inDName, dTab->dCount, dTab->dMinRes, dTab->dMaxRes );
           }



/***/
/* Adjust search limits.
/* Report status of input couplings.
/***/

         if (dcInfo->segFirst < -99999)   dcInfo->segFirst = dTab->dMinRes;
         if (dcInfo->segLast  < -99999)   dcInfo->segLast  = dTab->dMaxRes;

         if (dcInfo->segFirst > dcInfo->segLast) {ISWAP( dcInfo->segFirst, dcInfo->segLast );}

         if (dcInfo->aFirst < -99999)     dcInfo->aFirst = dcInfo->segFirst;


         (void) updatePDBAngles( pdbInfo );
         (void) updatePDBShifts( pdbInfo );
         (void) updatePDBSurf( pdbInfo );



/***/
/* Allocate and extract coord lists, etc:
/***/

         atomCount = pdbInfo->atomCount;
         atomList  = pdbInfo->atomList;


         xI    = dFltAlloc( "matrix", atomCount+1+MaxCountCharge );
         yI    = dFltAlloc( "matrix", atomCount+1+MaxCountCharge );
         zI    = dFltAlloc( "matrix", atomCount+1+MaxCountCharge );

         if (!xI || !yI || !zI)
            {
             FPR( stderr, "Error allocating coord buffers.\n" );
             error = 1;
             goto shutdown;
            }


         chargeAtomID    = intAlloc( "matrix", MaxCountCharge+1 );
         charges         = fltAlloc( "matrix", MaxCountCharge+1 );
         potential       = dFltAlloc( "matrix", MAX_DIG_POT+1 );
         ddots           = dFltAlloc( "matrix", (MAX_ENSEMBLE+1)*3 );

         if (!chargeAtomID || !charges || !potential || !ddots )
            {
             FPR( stderr, "Error allocating charge and potential buffers.\n" );
             error = 1;
             goto shutdown;
            }



/***/
/* Select the input shifts and couplings to fit:
/***/

         validCountD  = 0;
         validCountCS = 0;

         dObsList     = dTab->dObsList;
         dCalcList    = dTab->dCalcList;
         dCount       = dTab->dCount;

         validCountD  = mapDCInfo( pdbInfo, dTab );

         validCountD  = selectD( dObsList,
                            dCount,
                            dcInfo->segFirst,
                            dcInfo->segLast,
                            dcInfo->excludeMode );

         validCountD  = selTypDC( dObsList, dCount, atomNameS1, atomNameS2 );

/*
         FPR( stdout,
              "REMARK %d couplings selected (%s%s).\n",
              validCountD,
              dcInfo->excludeMode == DC_EXCLUDE     ? "exclusive" : "inclusive",
              dcInfo->excludeMode == DC_INCLUDE_MAX ? " max" : "" );
*/

         validCountD  = 
          buildDC( pdbInfo, dTab, dcInfo->aFirst-dcInfo->segFirst, dcInfo->segFirst, wList, wCount );

         FPR( stdout,
              "REMARK %d couplings selected (%s%s).\n",
              validCountD,
              dcInfo->excludeMode == DC_EXCLUDE     ? "exclusive" : "inclusive",
              dcInfo->excludeMode == DC_INCLUDE_MAX ? " max" : "" );



/***/
/* PDB Selection:
/***/

         if( sInfo.res1 == 0 && sInfo.resN == 0 )
           {
              sInfo.res1  = pdbInfo->minResID;
              sInfo.resN  = pdbInfo->maxResID;
           }

         selectCount = pdbSelect( pdbInfo, &sInfo );


         FPR( stderr, "REMARK PDB Selection:" );

         if (sInfo.selectMode == PDB_SELECT_ALL)
             FPR( stderr, " All Atoms.\n" );
         else
            FPR( stderr, 
                 " Residues %d to %d  Residues: %s  Atoms: %s  Chains %s  Segs: %s.\n",
                 sInfo.res1,
                 sInfo.resN,
                 sInfo.resNameList,
                 sInfo.atomNameList,
                 sInfo.chainNameList,
                 sInfo.segNameList );


         if( !palesInfo->Hflag )
           {
              sInfo.selectID     = 0;
              sInfo.atomNameMode = STRMATCH_GLOB;
              sInfo.selectMode   = PDB_SELECT_PARSE;
              (void) strcpy( sInfo.atomNameList, "[1-3]H* H*" );
 
              selectCount = pdbSelect( pdbInfo, &sInfo );
           }

 

/***/
/* Selection of mode of PALES:
/***/


         if ( palesInfo->palesMode > 6 )
           {

/***/
/***/
/* DC mode:
/***/
/***/

              if( !mapCount )
                {
                   if( error = simDC( pdbInfo, dTab, palesInfo, 
                                      remark, outDName, pdbRotName, 
                                      palesInfo->verbFlag ) )
                     {
                        FPR( stderr, "Matrix: error in best-fit of data.\n");
                        goto shutdown;
                     }
                }


              else
                {


/***/
/* 2 checks not done compared to orderten_svd.c 		*/
/* ----------------------------------------------------------	*/
/* A) 								*/
/*  for physicaly acceptable value the diagonal elements 	*/
/*  are between [-0.5-1] 					*/
/*  the offdiagonal elements are between [-0.75-0.75]		*/
/* B)								*/
/*  checking acceptibility of eigenvalues			*/
/***/

                   if( error = mapSimDC( pdbInfo, dTab, palesInfo, 
                                         outDName, pdbRotName,
                                         palesInfo->verbFlag, 
                                         outMapName, mapCount, remark ) )
                     {
                        FPR( stderr, "Matrix: error in best-fit of data.\n");
                        goto shutdown;
                     }

                } /*  end of error mapping  */

          }


/***/
/***/
/* PALES mode:
/***/
/***/

       else
          {


/***/
/* Read configuration file:
/* Read file with electrostatic potential:
/***/

            if (error = readConfig( inpName, palesInfo ) )
              {
                 FPR( stderr, "Error reading configuration file %s\n", inpName );
                 goto shutdown;
              }



/***/
/* Refine model of anisotropic medium
/***/

            refine_bicelle( palesInfo );



            switch( palesInfo->palesMode )
               {


                 case PALES_STERIC:

/***/
/* Setup orientation grid:
/***/

    density = palesInfo->densityNSC;
    if( UNSP(density,&dots,&n_dots) )
      {
        FPR( stderr, "Matrix: unsp error occured.\n");
	error = 1;
        goto shutdown;
      }
    for( i=0; i<n_dots; i++ )
       {
         ddots[i*3]   = dots[i*3];
         ddots[i*3+1] = dots[i*3+1];
         ddots[i*3+2] = dots[i*3+2];
       }


/***/
/* Select coordinates for simulation:
/***/

   if ( error = chCoorSelect( pdbInfo, palesInfo, &atomSelCount, xI, yI, zI, 
			chCount, chargeAtomID ) )
      {
        FPR( stderr, "pales: error reading selection of atoms.\n");
        goto shutdown;
      }

   FPR( stderr, "REMARK %d atoms selected for simulation.\n", atomSelCount);



/***/
/* Simulate macromolecular alignment:
/***/

   if ( error = simulateAlignment( palesInfo, atomSelCount, xI, yI, zI , 
		chCount, chargeAtomID, charges, n_dots, ddots,
		digPotential, potential, offset, stepsize, saupe ) )
       {
        FPR( stderr, "Matrix: Error simulating alignment.\n");
        goto shutdown;
       }



/***/
/* Selection of sign according to Frank's convention:
/***/

   for( i = 0; i < 5; i++ )
      {
         saupe[i] = (-1.0)*saupe[i];
      }



/***/
/* Analyze order matrix:
/***/

    if( analyzeSaupe(saupe, irrSaupe, &gMag, eigVec, eigVal, &da, &dr,
        alphaList, betaList, gammaList, palesInfo->eulerFlag) )
      {
        FPR( stderr, "Matrix: error in analyzeSaupe.\n");
        error = 1;
        goto shutdown;
       }



/***/
/* Analyze order matrix with different axes selection:
/***/

/*
    if( error = analyzeEigen( eigVal, eigVec, 0 ) )
      {
        FPR( stderr, "Matrix: error in analyzeEigen.\n");
        goto shutdown;
      }
*/



/***/
/* Compare predicted dipolar couplings with observed ones:
/* Analyze observed dipolar couplings:
/* Write output:
/***/

    dTab->dTensor.tensorMode = DC_STERIC;

    if( validCountD != 0 )
      {
        for( i = 0; i < 5; i++ )
          {
              dTab->dTensor.userSaupe[i] = saupe[i];
          }

         (void) fitD( pdbInfo, dTab, palesInfo->verbFlag );
      }
    else
     {
        for( i = 0; i < 5; i++ )
          {
              dTab->dTensor.saupe[i] = saupe[i];
          }
      }


    (void) writeDCmod( outDName, "w", dTab, palesInfo, 
                result );



/***/
/* Write molecule rotated such that alignment tensor coincide with laboratory frame:
/***/

   if( error = storePDBCoords( pdbInfo ) )
     {
        FPR( stderr, "matrix: error in storePDBCoords.\n" );
        goto shutdown;
     }


   rotAngle[0]=-1.*gammaList[0];
   rotAngle[1]=-1.*betaList[0];
   rotAngle[2]=-1.*alphaList[0];

   if( error = rotateDomain(rotAngle,pdbInfo,pdbInfo->minResID,pdbInfo->maxResID) )
     {
	FPR( stderr, "matrixt: error in rotateDomain." );
	goto shutdown;
     }



   if( *pdbRotName )
     {
        (void) writePDB( pdbRotName, pdbInfo, PDB_SELECT_NULL, remark );
     }


   if( error = restorePDBCoords( pdbInfo ) )
     {
         FPR( stderr, "matrix: error in storePDBCoords.\n" );
         goto shutdown;
     }


		   break;
 
              default:
		   break;
              } /* end of switch */


          } /* end of if palesInfo->palesMode <= 6 */


     } /* end of if palesInfo->palesMode < 8 */


/***/
/* Exit point:
/***/

shutdown:

    (void) pdbFree( pdbInfo );
    (void) dTabFree( dTab );
    (void) csTabFree( csTab );
    (void) fitAlignFree();
    (void) dipolarFree();
    (void) surfFree();

    (void) deAlloc( "listArg", wList, sizeof(float)*wCount );

    (void) showAlloc();


    (void) deAlloc( "matrix", pdbInfo,      sizeof( struct PDBInfo ));
    (void) deAlloc( "matrix", dTab,         sizeof( struct DTable ));
    (void) deAlloc( "matrix", palesInfo,    sizeof( struct PalesInfo ));
    (void) deAlloc( "matrix", dcInfo,       sizeof( struct DcInfo ));

    (void) deAlloc( "matrix", xI,      sizeof(double)*(atomCount+1+MaxCountCharge ) );
    (void) deAlloc( "matrix", yI,      sizeof(double)*(atomCount+1+MaxCountCharge ) );
    (void) deAlloc( "matrix", zI,      sizeof(double)*(atomCount+1+MaxCountCharge ) );

    (void) deAlloc( "matrix", chargeAtomID,      sizeof(int)*(MaxCountCharge+1) );
    (void) deAlloc( "matrix", charges,           sizeof(float)*(MaxCountCharge+1) );
    (void) deAlloc( "matrix", potential,         sizeof(double)*(MAX_DIG_POT+1) );
    (void) deAlloc( "matrix", ddots,             sizeof(double)*(MAX_ENSEMBLE+1)*3 );

    (void) dynDone();


    return( error );
}


/* =================================================================== */
/* =================================================================== */
/* =================================================================== */

/* =================================================================== */
int palesInfoNull( struct PalesInfo *palesInfo )
{
    (void) elInfoNull( &palesInfo->elInfo );
    (void) gridInfoNull( &palesInfo->gridInfo );

    palesInfo->palesMode       = PALES_STERIC;
    palesInfo->lcType          = 1;
    palesInfo->concentration   = 0.05;
    palesInfo->densityNSC      = 100;
    palesInfo->digPsiNSC       = 18;
    palesInfo->radiusModel     = 20.0;
    palesInfo->gridSpacing     = 0.2;
    palesInfo->lcOrder         = 0.8;     /* liquid crystall order: [0-1] */
    palesInfo->Hflag           = 0;
    palesInfo->surfFlag        = 1;       /* select surface accessible atoms [0/1] */
    palesInfo->radiusAtom      = 0.0;

    palesInfo->includeRim      = 0.0;
    palesInfo->molarRatio      = 3.0;
    palesInfo->lengthModel     = 200.0;
    palesInfo->orderDiscRim    = 0.25;

    palesInfo->rSingleCell     = 0.0;

    palesInfo->opCode          = OP_SUB;
    palesInfo->c1              = 1.0;
    palesInfo->c2              = 1.0;

    palesInfo->eulerFlag       = 1;
    palesInfo->verbFlag        = 1;


    return( 0 );
}
/* =================================================================== */

/* =================================================================== */
int dcInfoNull( struct DcInfo *dcInfo )
{
    dcInfo->segFirst    = -999999;
    dcInfo->segLast     = -999999;
    dcInfo->aFirst      = -999999;
    dcInfo->segID       = 1;
    dcInfo->segStep     = 0;

    dcInfo->excludeMode = DC_INCLUDE_MAX;

    return( 0 );
}
/* =================================================================== */

/* =================================================================== */
int gridInfoNull( struct GridInfo *gridInfo )
{
    gridInfo->densE           = 200;      /* density of dots on ensemble sphere */
    gridInfo->angleDomFlag    = 0;        /* change rel. orientation of domains (0/1) */
    gridInfo->angleDom[0]     = 0.0;      /* Euler angles for adjusting */
                                          /* angle between 2 domains */
    gridInfo->angleDom[1]     = 0.0; 	  
    gridInfo->angleDom[2]     = 0.0; 	 
    gridInfo->resDom1Origin   = 0;        /* origin for orientation of domain 1 */
    gridInfo->resDom2Origin   = 0;        /* origin for orientation of domain 2 */
    gridInfo->startResDom2    = 0;        /* domain boundaries */
    gridInfo->startResDom1    = 0;
    gridInfo->endResDom1      = 0;
    gridInfo->endResDom2      = 0;
    gridInfo->digPsiE         = 1;        /* digitization of rotation about */
				          /* 3rd Euler angle */
    gridInfo->lBetaCone       = 0.0;
    gridInfo->hBetaCone       = 30.0;
    gridInfo->cAlphaCone      = 0.0;
    gridInfo->errAlphaCone    = 180.0;
    gridInfo->startResPart1   = 0;
    gridInfo->startResPart2   = 0;
    gridInfo->endResPart1     = 0;
    gridInfo->endResPart2     = 0;
    gridInfo->clashFlag       = 1;        /* check steric clash of domains (0/1) */
    gridInfo->clashDist       = 3.6;      /* distance for steric clash */
    gridInfo->initOrientFlag  = 1;
    gridInfo->pivotCenter     = 1;        /* 1: translate molecule to center of mass */
				          /* prior to PALES prediction */
    gridInfo->fixRotDomain    = 0;        /* 1: rotate protein such that domain which */
				          /* is rotated is returned to starting position */

    return( 0 );
}
/* =================================================================== */

/* =================================================================== */
int elInfoNull( struct ElInfo *elInfo )
{
    elInfo->elFlag          = 0;
    elInfo->pH              = 7.0;
    elInfo->concNaCl        = 0.2;        /* [M] molarity mol/l */
    elInfo->pkaFlag         = 0;          /* using pKa default values */
    elInfo->chTerminiFlag   = 1;          /* flag for putting charges at termini (0/1)*/
    elInfo->chDensSurf      = -0.47;      /* surface charge density [e/nm2] */
    elInfo->kEnergy         = 1.0;        /* scaling for EL energy */
    elInfo->kPot            = 1.0;        /* scaling for EL potential */
    elInfo->cutoffPot       = 400.0;      /* el cutoff [A] */
    elInfo->kNegCh          = 1.0;
    elInfo->kPosCh          = 1.0;
    elInfo->temperaturePot  = 298.150;    /* temperature [K] */
    elInfo->chFlag          = 0;

    return( 0 );
}
/* =================================================================== */

/* ======================== rotate vector ================================= */
void rotate_vector(double theta, double phi, double psi, double *input_vector,
        double *rotated_vector)
/* rotation of a vector according to specified Euler angles */
{
        int i,j;
        double a,b,c;
        double Dreh_M[4][4];


/*
        for(i=1;i<=3;i++)
                rotated_vector[i] = 0;
        rotated_vector[1]=( C(phi)*C(psi)-S(phi)*C(theta)*S(psi) )*input_vector[1];
        rotated_vector[1]+=( S(phi)*C(psi)+C(phi)*C(theta)*S(psi) )*input_vector[2];
        rotated_vector[1]+=( S(theta)*S(psi) )*input_vector[3];

        rotated_vector[2]=( -1*C(phi)*S(psi)-S(phi)*C(theta)*C(psi) )*input_vector[1];
        rotated_vector[2]+=( -1.*S(phi)*S(psi)+C(phi)*C(theta)*C(psi) )*input_vector[2];
        rotated_vector[2]+=( S(theta)*C(psi) )*input_vector[3];

        rotated_vector[3]=( S(phi)*S(theta) )*input_vector[1];
        rotated_vector[3]+=( -1*C(phi)*S(theta) )*input_vector[2];
        rotated_vector[3]+=C(theta)*input_vector[3];
*/

        /* rotation matrix as obtained from Matlab program */
        a=phi;
        b=theta;
        c=psi;
        Dreh_M[1][1]=cos(a)*cos(b)*cos(c)-sin(a)*sin(c);
        Dreh_M[1][2]=sin(a)*cos(b)*cos(c)+cos(a)*sin(c);
        Dreh_M[1][3]=-sin(b)*cos(c);
        Dreh_M[2][1]=-cos(a)*cos(b)*sin(c)-sin(a)*cos(c);
        Dreh_M[2][2]=-sin(a)*cos(b)*sin(c)+cos(a)*cos(c);
        Dreh_M[2][3]=sin(b)*sin(c);
        Dreh_M[3][1]=cos(a)*sin(b);
        Dreh_M[3][2]=sin(a)*sin(b);
        Dreh_M[3][3]=cos(b);

        /* rotation matrix as defined by Korn & Korn 1961 pg. 475 */
/*
        Dreh_M[1][2]=-1*Dreh_M[1][2];
        Dreh_M[2][1]=-1*Dreh_M[2][1];
        Dreh_M[3][1]=-1*Dreh_M[3][1];
*/
        for(i=1;i<=3;i++) {
                rotated_vector[i] = 0;
                for(j=1;j<=3;j++) rotated_vector[i] += Dreh_M[i][j]*input_vector[j];
        }
}
/* ========================================================================= */


/* ================== rotate coordinate system =========================== */
void rotate_coorsystem(double theta, double phi, double psi, double *x_axis_tensor,
        double *y_axis_tensor, double *z_axis_tensor)
/* Rotation of axis of alignment tensor according to specified Euler angles */
{
        double x_axis_xray[4] = {0.0,1.0,0.0,0.0};
        double y_axis_xray[4] = {0.0,0.0,1.0,0.0};
        double z_axis_xray[4] = {0.0,0.0,0.0,1.0};

        rotate_vector(theta, phi, psi, x_axis_xray, x_axis_tensor);
        rotate_vector(theta, phi, psi, y_axis_xray, y_axis_tensor);
        rotate_vector(theta, phi, psi, z_axis_xray, z_axis_tensor);
}
/* ========================================================================= */


/* ============ simulate sterically induced alignment ======================= */
void ssia( struct PalesInfo *palesInfo, int no_atoms, 
           double *x, double *y, double *z,
           double max_protein_length, double min_protein_length,
           int n_dots, double *dots, double *saupe )
{
  int     i, j, k, pts, add_pts;
  int     total_bic=0, total_ph=0, total_orient;
  int     no_ok_pts_bic=0, no_ok_pts_ph=0;
  double  atom_vec[4];
  double  radius=0.0,radius2=0.0;
  double  scale_fac=0.0;
  static  double  saupe_dis[MAX_NO_PT][5];
  static  double  saupe_bic[5], saupe_ph[5];
  int     m=0, orientation_flag=1;

  double  x_axis_tensor[4],y_axis_tensor[4],z_axis_tensor[4];
  double  Dreh_m[4][4];
  float   *x_temp,*y_temp,*pot_dist,*grid_el;
  double  *total;

  int     a,lengthGrid;
  int     t_pts, no_pt_ph;
  int     no_add_orient, no_add_ph_orient;
  float   offset_grid=0;
  double  d,p;
  double  end_dist, density=1.46;


  int     kk,mm;
  double  inc_psi;
  double  beta,gamma,alpha;

  int     o_pts;
  double  szz,sxx_yy,sxy,sxz,syz;
  double  tensor_xaxis[4],tensor_yaxis[4],tensor_zaxis[4];

  double  saupe_dev[5];
  double  sum1=0.0,sum2=0.0,sum3=0.0,sum4=0.0,sum5=0.0,sum6=0.0;
  FILE    *fp0,*fp1,*fp2,*fp3,*fp4;


  x_temp   = (float  *)NULL;
  y_temp   = (float  *)NULL;
  pot_dist = (float  *)NULL;
  grid_el  = (float  *)NULL;
  total    = (double *)NULL;


  x_temp = (float *)malloc((size_t) (no_atoms+1)*sizeof(float));
        if (!x_temp) nrerror("allocation failure for x_temp");
  y_temp = (float *)malloc((size_t) (no_atoms+1)*sizeof(float));
        if (!y_temp) nrerror("allocation failure for y_temp");
  pot_dist = (float *)malloc((size_t) (MAX_NO_PT+1)*sizeof(float));
        if (!pot_dist) nrerror("allocation failure for pot_dist");
  grid_el = (float *)malloc((size_t) (MAX_NO_PT)*sizeof(float));
        if (!grid_el) nrerror("allocation failure for grid_el");
  total = (double *)malloc((size_t) 2*(MAX_NO_PT)*sizeof(double));
        if (!total) nrerror("allocation failure for total");
 


  radius = palesInfo->radiusAtom;


/****************/
/* theta = beta */
/* phi = alpha  */
/* psi = gamma  */
/****************/

/* the Euler angles theta and phi are equivalent to the polar angles */
/* which determine the orientation of the z-axis */


  /**************************************************/
  /*** check uniformity of orientational sampling ***/
  inc_psi=360.0/palesInfo->digPsiNSC;
/*
  fp0 = fopen("zaxis.txt", "w");
  fp1 = fopen("beta.txt", "w");
  fp2 = fopen("alpha.txt", "w");
  fp3 = fopen("gamma.txt", "w");
  fp4 = fopen("xaxis.txt", "w");
  printf("start test\n");
*/
  for (i=0; i<n_dots; i++){

        /* mapping from three-dimensional Cartesian coordinates */
        /* to spherical coordinates as in cart2sph.m */
        alpha = atan2(dots[i*3+1],dots[i*3]);
        beta = atan2(dots[i*3+2],sqrt(dots[i*3]*dots[i*3]+dots[i*3+1]*dots[i*3+1]));
        beta = PI/2.0-beta;

        for(mm=0;mm<palesInfo->digPsiNSC;mm++) {                   
                gamma = RAD(0.0+mm*inc_psi);

                rotate_coorsystem(beta,alpha,gamma,tensor_xaxis,
                        tensor_yaxis,tensor_zaxis);

/*
                fprintf(fp0,"%.3lf %.3lf %.3lf\n",
                        tensor_zaxis[1],tensor_zaxis[2],tensor_zaxis[3]);
                fprintf(fp4,"%.3lf %.3lf %.3lf\n",
                        tensor_xaxis[1],tensor_xaxis[2],tensor_xaxis[3]);
                fprintf(fp1,"%.3lf %s\n",DEG(beta),"1");
                fprintf(fp2,"%.3lf %s\n",DEG(alpha),"1");
                fprintf(fp3,"%.3lf %s\n",DEG(gamma),"1");

                beta  = tensor_zaxis[1]*tensor_xaxis[1]+tensor_zaxis[2]*tensor_xaxis[2]+tensor_zaxis[3]*tensor_xaxis[3];
                alpha = tensor_zaxis[1]*tensor_yaxis[1]+tensor_zaxis[2]*tensor_yaxis[2]+tensor_zaxis[3]*tensor_yaxis[3];
                gamma = tensor_xaxis[1]*tensor_yaxis[1]+tensor_xaxis[2]*tensor_yaxis[2]+tensor_xaxis[3]*tensor_yaxis[3];
                if( beta>10e-16 || alpha>10e-16 || gamma>10e-16 )
                        printf("\nnorm: %e   %e   %e\n\n",beta,alpha,gamma);
*/

                /* calculate elements of Saupe matrix for Euler angles i, */
                /* e.g. polar angles theta (beta) and phi (alpha) */
                sum1 += Szz(tensor_zaxis);
                sum2 += Sxx(tensor_xaxis);
                sum3 += Syy(tensor_yaxis);
                sum4 += Sxy(tensor_xaxis,tensor_yaxis);
                sum5 += Sxz(tensor_xaxis,tensor_zaxis);
                sum6 += Syz(tensor_yaxis,tensor_zaxis);
        } /* end gamma loop */
  } /* end surface dot loop */
/*
  printf("end test\n");
  fclose(fp0);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
*/

  o_pts = n_dots*palesInfo->digPsiNSC;
/*
  saupe_dev[0]=sum1/(1.0*o_pts);
  saupe_dev[1]=(sum2-sum3)/(1.0*o_pts);
  saupe_dev[2]=sum4/(1.0*o_pts);
  saupe_dev[3]=sum5/(1.0*o_pts);
  saupe_dev[4]=sum6/(1.0*o_pts);
  for (i=0; i<5; i++)
        printf("o_pts: %d    saupe_dev[%d]: %e\n",o_pts,i,saupe_dev[i]);
*/



  /*************/
  /*** setup initial tranlational grid ***/
   for( i=0; i<MAX_NO_PT; i++ )   grid_el[i] = 0.0;

   offset_grid = min_protein_length + palesInfo->radiusModel;
   t_pts       = (max_protein_length-min_protein_length)/palesInfo->gridSpacing;
   if( t_pts >= MAX_NO_PT )
     {
       t_pts =  MAX_NO_PT - 1;
       FPR( stderr, "WARNING: concentration too low - magnitude of alignment tensor incorrect\n" );
     }
   for( i=0; i<=t_pts; i++ )
        grid_el[i] = offset_grid + i * palesInfo->gridSpacing;
   no_pt_ph = t_pts;



  /*************/
  /*** calculate scaling ***/
   if( palesInfo->lcType == 1 )
        no_add_orient=(palesInfo->rSingleCell-grid_el[t_pts])/palesInfo->gridSpacing;
   if( palesInfo->lcType == 2 )
     {
        end_dist = palesInfo->radiusModel/sqrt(palesInfo->concentration/density); 
        a        = (end_dist-palesInfo->radiusModel) / palesInfo->gridSpacing;
        no_add_ph_orient = a + (palesInfo->gridSpacing/palesInfo->radiusModel) * (a*(a+1)/2.0);
     }



  /*************/
  /*** plane ***/
   if(palesInfo->lcType==1) {

        /* printf("Starting with bicelle simulation ...\n"); */


        /*** initialization ***/
        for( j=0; j<=t_pts; j++ )   pot_dist[j] = radius-(grid_el[j]-palesInfo->radiusModel);

        for( i=0; i<5; i++)
                for( j=0; j<=t_pts; j++ )
                        saupe_dis[j][i] = 0.0;
        for( j=0; j<=t_pts; j++ )
                total[j]=0;



        /*** perform simulation ***/
        kk=0;
        for(i=0; i<n_dots; i++){

           /* mapping from three-dimensional Cartesian coordinates */
           /* to spherical coordinates as in cart2sph.m */
           alpha = atan2(dots[i*3+1],dots[i*3]);
           beta = atan2(dots[i*3+2],sqrt(dots[i*3]*dots[i*3]+dots[i*3+1]*dots[i*3+1]));
           beta = PI/2.0-beta;

           for(mm=0;mm<palesInfo->digPsiNSC;mm++) {                        
                gamma = RAD(0.0+mm*inc_psi);

                /* report how many orientations have been tested */
                /* if(fmod(++kk,1000)==0.0)     printf("# of orientations: %d\n",kk); */
                /* if(fmod(++kk,200)==0.0)      printf(".\n"); */


                /* calculate elements of Saupe matrix for Euler angles i, */
                /* e.g. polar angles theta (beta) and phi (alpha) */
                rotate_coorsystem(beta,alpha,gamma,tensor_xaxis,
                        tensor_yaxis,tensor_zaxis);

                szz = Szz(tensor_zaxis);
                sxx_yy = Sxx(tensor_xaxis)-Syy(tensor_yaxis);
                sxy = Sxy(tensor_xaxis,tensor_yaxis);
                sxz = Sxz(tensor_xaxis,tensor_zaxis);
                syz = Syz(tensor_yaxis,tensor_zaxis);


                /* rotate coordinate system (inverse rotation) */
                rotate_coorsystem(-1.*beta,-1.*gamma,-1.*alpha,x_axis_tensor,
                        y_axis_tensor,z_axis_tensor);
                for(k=1;k<=3;k++){
                        Dreh_m[k][1]=x_axis_tensor[k];
                        Dreh_m[k][2]=y_axis_tensor[k];
                        Dreh_m[k][3]=z_axis_tensor[k];
                }


                /* calculate coordinates in rotated frame */
                for(k=1;k<=no_atoms;k++) {
                        atom_vec[1]=x[k];       /* select cartesian */
                        atom_vec[2]=y[k];       /* coordinates of */
                        atom_vec[3]=z[k];       /* atom k */

                        x_temp[k]=0;
                        for(j=1;j<=3;j++)       x_temp[k]+=Dreh_m[j][1]*atom_vec[j];
                }
                

                /* sort x coordinates into ascending order */
                /* slows the program down */
                /* shell(no_atoms,x_temp); */


                /* check hard potential for different distances */
                for(j=0;j<=t_pts;j++) {

                        orientation_flag=1;
                        for(m=1;m<=no_atoms;m++) {
                                if(x_temp[m] <= pot_dist[j]) { 
                                        orientation_flag=0;
                                        break;
                                }
                        } 

                        if(orientation_flag) {
                                saupe_dis[j][0] += szz;
                                saupe_dis[j][1] += sxx_yy;
                                saupe_dis[j][2] += sxy;
                                saupe_dis[j][3] += sxz;
                                saupe_dis[j][4] += syz;
                                total[j]++;
                        }       
                }

           } /* end gamma loop */
        } /* end surface dot loop */


        /* refine and norm Saupe matrix for specific translational position */
        for(i=0;i<5;i++)
                for(j=0;j<=t_pts;j++) {
                        if(total[j]>0)
                                saupe_dis[j][i] = saupe_dis[j][i] / (1.0*total[j]) - saupe_dev[i];
                        else
                                saupe_dis[j][i] = 0.0;
                }


        no_ok_pts_bic=0;
        total_bic=0;
        for(j=0;j<=t_pts;j++) {
                if(total[j]!=0) no_ok_pts_bic++;
                total_bic += total[j];
        }


        /* calculate Saupe matrix for range of grid points */
        /* over which simulation is performed */
        for(i=0;i<5;i++) {
                saupe_bic[i] = 0;
                for(j=0;j<=t_pts;j++) {
                        saupe_bic[i] += (total[j]*1.0/(o_pts*1.0)) * saupe_dis[j][i]; 
                }
        }

  }  /* end bicelle part */




  /****************/
  /*** cylinder ***/
   if(((palesInfo->lcType==1) && (palesInfo->includeRim==1.0)) || (palesInfo->lcType==2)) {
        
       /* printf("Starting with phage simulation ...\n"); */


        /* initialization */
        radius2 = palesInfo->radiusModel + radius;


        /* set up complete grid */
        end_dist = palesInfo->radiusModel/sqrt(palesInfo->concentration/density)-palesInfo->radiusModel; 
        i=0;
        d=grid_el[no_pt_ph];
        while( d <= end_dist )
             {  
                i++;
                if( no_pt_ph+i == MAX_NO_PT )
                  {
                      FPR( stderr, "WARNING: concentration too low - magnitude of alignment tensor incorrect\n" );
                      break;
                  }

                grid_el[no_pt_ph+i] = grid_el[no_pt_ph] + palesInfo->gridSpacing*i;
                d= grid_el[no_pt_ph+i];
             }
        lengthGrid = i+no_pt_ph;
        /* lengthGrid = no_pt_ph; */


        for(i=0;i<5;i++)
                for(j=0;j<=lengthGrid;j++)
                        saupe_dis[j][i] = 0.0;
        for(j=0;j<=lengthGrid;j++)
                total[j]=0.0;



        /*** perform simulation ***/
        kk=0;
        for(i=0; i<n_dots; i++){

           /* mapping from three-dimensional Cartesian coordinates */
           /* to spherical coordinates as in cart2sph.m */
           alpha = atan2(dots[i*3+1],dots[i*3]);
           beta = atan2(dots[i*3+2],sqrt(dots[i*3]*dots[i*3]+dots[i*3+1]*dots[i*3+1]));
           beta = PI/2.0-beta;

           for(mm=0;mm<palesInfo->digPsiNSC;mm++) {                        
                gamma = RAD(0.0+mm*inc_psi);

                /* report how many orientations have been tested */
                /* if(fmod(++kk,1000)==0.0)     printf("# of orientations: %d\n",kk); */
                /* if(fmod(++kk,200)==0.0)      printf(".\n");  */


                /* calculate elements of Saupe matrix for Euler angles i, */
                /* e.g. polar angles theta (beta) and phi (alpha) */
                rotate_coorsystem(beta,alpha,gamma,tensor_xaxis,
                        tensor_yaxis,tensor_zaxis);

                szz = Szz(tensor_zaxis);
                sxx_yy = Sxx(tensor_xaxis)-Syy(tensor_yaxis);
                sxy = Sxy(tensor_xaxis,tensor_yaxis);
                sxz = Sxz(tensor_xaxis,tensor_zaxis);
                syz = Syz(tensor_yaxis,tensor_zaxis);


                /* rotate coordinate system (inverse rotation) */
                rotate_coorsystem(-1.*beta,-1.*gamma,-1.*alpha,x_axis_tensor,
                        y_axis_tensor,z_axis_tensor);
                for(k=1;k<=3;k++){
                        Dreh_m[k][1]=x_axis_tensor[k];
                        Dreh_m[k][2]=y_axis_tensor[k];
                        Dreh_m[k][3]=z_axis_tensor[k];
                }


                /* calculate coordinates in rotated frame */
                for(k=1;k<=no_atoms;k++) {
                        atom_vec[1]=x[k];       /* select cartesian */
                        atom_vec[2]=y[k];       /* coordinates of */
                        atom_vec[3]=z[k];       /* atom k */

                        x_temp[k]=0;
                        y_temp[k]=0;
                        for(j=1;j<=3;j++){
                                x_temp[k]+=Dreh_m[j][1]*atom_vec[j];
                                y_temp[k]+=Dreh_m[j][2]*atom_vec[j];
                        }
                }


                /* shift center of mass of protein to new position on grid */
                for(j=0;j<=lengthGrid;j++) {

                        p=1.0;
                        orientation_flag=1;
                        if(j<=no_pt_ph)
                        for(m=1;m<=no_atoms;m++) {
                                if( sqrt((x_temp[m]+grid_el[j])*(x_temp[m]+grid_el[j])+(y_temp[m]*y_temp[m])) <= radius2) { 
                                        orientation_flag=0;
                                        p=0.0;
                                        break;
                                }
                        }

                        saupe_dis[j][0] += p*szz;
                        saupe_dis[j][1] += p*sxx_yy;
                        saupe_dis[j][2] += p*sxy;
                        saupe_dis[j][3] += p*sxz;
                        saupe_dis[j][4] += p*syz;

                        total[j] += p;  
                }


           } /* end gamma loop */
        } /* end dots loop */



        /* refine and norm Saupe matrix for specific translational position */
        for(i=0;i<5;i++)
                for(j=0;j<=lengthGrid;j++) {
                        if(total[j]>0.0)
                                saupe_dis[j][i] = saupe_dis[j][i] / total[j] - saupe_dev[i];
                        else
                                saupe_dis[j][i] = 0.0;
                }



        no_ok_pts_ph=0;
        total_ph=0;
        for(j=0;j<=lengthGrid;j++) {
                if(total[j]!=0.0)       no_ok_pts_ph++;
                total_ph += total[j];
        }

        /* for(i=0;i<=lengthGrid;i++) printf("grid_el[%d]: %.1f  total[%d]: %.1f     saupe_dis[%d][0]: %e\n",
                                        i,grid_el[i],i,total[i],i,(1.0/(o_pts*1.0)) *saupe_dis[i][0]); */

        /* calculate Saupe matrix for range of grid points */
        /* for which simulation is performed */
        for(i=0;i<5;i++) {
                saupe_ph[i] = 0;
                for(j=0;j<=lengthGrid;j++) {
                        d = (total[j]/(o_pts*1.0)) * saupe_dis[j][i] * grid_el[j]/palesInfo->radiusModel;
                        /* if(i==0) printf("grid_el[%d]/palesInfo->radiusModel: %.1f   saupe_dis[%d][0]: %e   saupe_disSCALE[%d][0]: %e\n",
                                j,grid_el[j]/palesInfo->radiusModel,j,saupe_dis[j][0],j,d); */
                        saupe_ph[i] += d;
                }
        }

    } /* end phage part */


/*
    printf("Saupe bicelle\n");
    for(i=0;i<5;i++) printf("%e  ",saupe_bic[i]);
    printf("\n");
    printf("Saupe phage\n");
    for(i=0;i<5;i++) printf("%e  ",saupe_ph[i]);
    printf("\n");
*/

    /* calculate Saupe matrix according to employed model system */
    if ((palesInfo->lcType==1) && (palesInfo->includeRim==1.0))
        for(i=0;i<5;i++) {
                saupe[i] = 1.0/(1+(1.0/palesInfo->molarRatio)) * (saupe_bic[i] -  
                (1.0/palesInfo->molarRatio)*((1.0*no_ok_pts_bic) / (1.0*no_ok_pts_ph)) * palesInfo->orderDiscRim * saupe_ph[i]);
        }
    if ((palesInfo->lcType==1) && (palesInfo->includeRim!=1.0))
        for(i=0;i<5;i++) {
                saupe[i] = saupe_bic[i];
                /* saupe[i] = 1.0/(1+(1.0/palesInfo->molarRatio)) * (1-(1.0/palesInfo->molarRatio) * palesInfo->orderDiscRim) * saupe_bic[i];
*/
        }
    if (palesInfo->lcType==2)
        for(i=0;i<5;i++) {
                saupe[i] = saupe_ph[i];
        }


    if (palesInfo->lcType==2)
      {
        add_pts = no_add_ph_orient;
      }
    else
     {
        add_pts = total_bic/o_pts + no_add_orient;
     }


    for(i=0;i<5;i++)
      {
        saupe[i] = saupe[i] / (add_pts*scale_fac);
      }


  free(x_temp);
  free(y_temp);
  free(pot_dist);
  free(grid_el);
  free(total);
}
/* ========================================================================= */


/* ========================================================================= */
int coorSelect( struct PDBInfo *pdbInfo, struct PalesInfo *palesInfo,
                int *atomSelCount, 
                double *xI, double *yI, double *zI )
{
  int i,k,error;


/***/
/* Initialization:
/***/

  error = 0;



/***/
/* Select surface accessible atoms:
/***/

  if( palesInfo->surfFlag )
    {
      if( error = callNSC( pdbInfo ) )
        {
          goto shutdown;
        }
    }



/***/
/* Store selected coordinates in 1D arrays
/***/

  k = 0; 
  for( i = 0; i < pdbInfo->atomCount; i++ )
       {
       if( pdbInfo->atomList[i].selected )
          {
           k++;
           xI[k] = pdbInfo->atomList[i].x; 
           yI[k] = pdbInfo->atomList[i].y; 
           zI[k] = pdbInfo->atomList[i].z; 
          }
       }
  *atomSelCount = k;



/***/
/* Exit point:
/***/

shutdown:

     return( error );
}
/* ========================================================================= */


/* ========================================================================= */
int refine_bicelle( struct PalesInfo *palesInfo )
{
    int error;
        double mol_DHPC_free = 5;
        double mw_DMPC=677.94, mw_DHPC=453.51;
        double mol_DHPC,mol_DMPC;
        double wv_DHPC,wv_DMPC;
        double concentration_DMPC,concentration_DHPC;
        double molar_ratio_eff,volume_bic,no_bic_mil;
        double cube_bic, height_cylinder, wv_bic;
        double LZahl=6.022*1e+23;



    error = 0;


        /*** w/v concentration of DMPC and of DHPC ***/
        mol_DHPC = palesInfo->concentration * 1000 / (palesInfo->molarRatio*mw_DMPC+mw_DHPC);
        mol_DMPC = palesInfo->molarRatio*mol_DHPC;
        wv_DHPC = mol_DHPC*mw_DHPC;
        wv_DMPC = mol_DMPC*mw_DMPC;
        concentration_DMPC = wv_DMPC/1000;
        concentration_DHPC = wv_DHPC/1000;


        /*** effective molar ratio DMPC:DHPC due to free DHPC ***/
        molar_ratio_eff = mol_DMPC / (mol_DHPC-mol_DHPC_free/1000);


        /*** density of bicelles taking into account free DHPC ***/
        if(palesInfo->concentration<0.01)  mol_DHPC_free=0;
        wv_bic = (mol_DMPC * mw_DMPC + (mol_DHPC-mol_DHPC_free/1000) * mw_DHPC) / 1000;


        /*** radius of bicelle employing effective molar ratio in bicelle ***/
        palesInfo->lengthModel = 0.5*palesInfo->radiusModel*molar_ratio_eff*(PI+sqrt(PI*PI+8/molar_ratio_eff));


        /*** distance to next bicelle** */
        /* volume of a single bicelle */
        volume_bic = (PI*palesInfo->lengthModel*palesInfo->lengthModel*2*palesInfo->radiusModel+0.5*2*PI*PI*palesInfo->lengthModel*palesInfo->radiusModel*palesInfo->radiusModel)*1e-24;
        /* # of bicelles/ml */
        /* assuming density of bicelle is 1, i.e. w/v is same as v/v */
        /* true bicelle density more like 1.15 */
        no_bic_mil = (1/volume_bic)*wv_bic;     /* # of bicelles/ml */
        /* no_bic_mil = (1/volume_bic)*palesInfo->concentration; */

        cube_bic = (1/no_bic_mil)*1e+24;        /* volume associated with one bicelle */
        /* model: bicelle in cylinder with radius (palesInfo->lengthModel+palesInfo->radiusModel) */
        /* --> height of cylinder */


        /* palesInfo->rSingleCell = palesInfo->radiusModel/palesInfo->concentration; */
        /* palesInfo->rSingleCell = palesInfo->radiusModel/concentration_DMPC; */
        palesInfo->rSingleCell = palesInfo->radiusModel/wv_bic;


/*
        palesInfo->rSingleCell *= (1-mol_DHPC_free/1000*mw_DHPC/1000);
        palesInfo->rSingleCell *= (1-mol_DHPC_free/1000*LZahl/1000 * (2./3.)*palesInfo->radiusModel*PI*(5*5)*1e-24);
        printf("palesInfo->rSingleCell: %e\n",palesInfo->rSingleCell);
*/

        /*** order for rim, relative to order for center ***/
        palesInfo->orderDiscRim=PI*palesInfo->lengthModel/(4*(PI*palesInfo->lengthModel+2*palesInfo->radiusModel));

/*
        printf("mol_DHPC_free/1000*LZahl/1000 * (2./3.)*palesInfo->radiusModel*PI*(5*5)*1e-24 = %e\n",
                mol_DHPC_free/1000/1000*LZahl * (2./3.)*palesInfo->radiusModel*PI*(5*5)*1e-24);
        printf("100-mol_DHPC_free*mw_DHPC/1000 = %e\n",100-mol_DHPC_free*mw_DHPC/1000);
        printf("palesInfo->rSingleCell: %e\n",palesInfo->rSingleCell);
        printf("mol_DMPC: %.1f [mMol]  mol_DHPC: %.1f [mMol]\n",1000*mol_DMPC,1000*mol_DHPC);
        printf("DMPC: %.1f [mg/ml]  DHPC: %.1f [mg/ml]\n",wv_DMPC,wv_DHPC);
        printf("wv_bic: %e\n",wv_bic);
        printf("volume_bic: %e\n",volume_bic);
        printf("no_bic_mil: %e\n",no_bic_mil);
        printf("cube_bic: %e\n",cube_bic);
        printf("height_cylinder: %e\n",height_cylinder);
        printf("palesInfo->rSingleCell: %e\n",palesInfo->rSingleCell);
        printf("effective molar ratio in bicelle: %.1lf\n",molar_ratio_eff);
        printf("effective concentration of DMPC: %.1lf [percentage]\n",concentration_DMPC*100);
        printf("effective concentration of DHPC: %.1lf [percentage]\n",concentration_DHPC*100);
        printf("radius of ideal bicelle: %.1lf [A]\n",palesInfo->lengthModel);
        printf("order parameter of rim relative to center: %.2lf\n",palesInfo->orderDiscRim);
*/

        /*** return values ***/ 
        /* palesInfo->concentration=concentration_DMPC; */
        palesInfo->molarRatio = molar_ratio_eff;

    return error;
}
/* ========================================================================= */


/* ========================================================================= */
int make_eigVecList(double eigVec[3][3], double eigVecXlist[NO_LIST][3],
        double eigVecYlist[NO_LIST][3], double eigVecZlist[NO_LIST][3])
{
   int i,j,k=0;
   double eigVec_temp[3][3], vector_product[3];


   /*** initialize ***/
   for(i=0;i<4;i++)
        for(j=0;j<3;j++){
                eigVecXlist[i][j]=0.0;
                eigVecYlist[i][j]=0.0;
                eigVecZlist[i][j]=0.0;
        }


   /*** original eigenvectors ***/
   for(i=0;i<3;i++) {
        eigVec_temp[i][0]=eigVec[i][0];
        eigVec_temp[i][1]=eigVec[i][1];
        eigVec_temp[i][2]=eigVec[i][2];
   }
   cross_productM2(1,2,eigVec_temp,vector_product);
   for(i=0;i<3;i++)     eigVec_temp[i][0]=vector_product[i];
   k += mirror_solutions(k,eigVec_temp,eigVecXlist,eigVecYlist,eigVecZlist);


   /*** (-y x z)/-y/z ***/
   for(i=0;i<3;i++) {
        eigVec_temp[i][0]=eigVec[i][0];
        eigVec_temp[i][1]=-1*eigVec[i][1];
        eigVec_temp[i][2]=eigVec[i][2];
   }
   cross_productM2(1,2,eigVec_temp,vector_product);
   for(i=0;i<3;i++)     eigVec_temp[i][0]=vector_product[i];
   k += mirror_solutions(k,eigVec_temp,eigVecXlist,eigVecYlist,eigVecZlist);


   /*** (y x -z)/y/-z ***/
   for(i=0;i<3;i++) {
        eigVec_temp[i][0]=eigVec[i][0];
        eigVec_temp[i][1]=eigVec[i][1];
        eigVec_temp[i][2]=-1*eigVec[i][2];
   }
   cross_productM2(1,2,eigVec_temp,vector_product);
   for(i=0;i<3;i++)     eigVec_temp[i][0]=vector_product[i];
   k += mirror_solutions(k,eigVec_temp,eigVecXlist,eigVecYlist,eigVecZlist);


   /*** (-y x -z)/-y/-z ***/
   for(i=0;i<3;i++) {
        eigVec_temp[i][0]=eigVec[i][0];
        eigVec_temp[i][1]=-1*eigVec[i][1];
        eigVec_temp[i][2]=-1*eigVec[i][2];
   }
   cross_productM2(1,2,eigVec_temp,vector_product);
   for(i=0;i<3;i++)     eigVec_temp[i][0]=vector_product[i];
   k += mirror_solutions(k,eigVec_temp,eigVecXlist,eigVecYlist,eigVecZlist);

  return( k );
}
/* ========================================================================= */


/* =============================================================== */
int cross_productM2(int n, int m, double eigVec[3][3], double product[3])
{
   product[0]=eigVec[1][n]*eigVec[2][m]-eigVec[2][n]*eigVec[1][m];
   product[1]=eigVec[2][n]*eigVec[0][m]-eigVec[0][n]*eigVec[2][m];
   product[2]=eigVec[0][n]*eigVec[1][m]-eigVec[1][n]*eigVec[0][m];

   return( 1 );
}
/* =============================================================== */


/* =============================================================== */
int add_eigVec_List(double eigVecXlist[NO_LIST][3], double eigVecYlist[NO_LIST][3],
        double eigVecZlist[NO_LIST][3], int n, double eigVec[3][3])
{
   int i;


   for(i=0;i<3;i++) {
        eigVecXlist[n][i]=eigVec[i][0];
        eigVecYlist[n][i]=eigVec[i][1];
        eigVecZlist[n][i]=eigVec[i][2];
   }

   return( 1 );
}
/* =============================================================== */


/* =============================================================== */
int mirror_solutions(int n, double eigVec[3][3], double eigVecXlist[NO_LIST][3],
        double eigVecYlist[NO_LIST][3], double eigVecZlist[NO_LIST][3])
{
   int i;
   double eigVec_temp[3][3];


   /*** initialize ***/
   for(i=0;i<3;i++) {
        eigVec_temp[i][0]=eigVec[i][0];
        eigVec_temp[i][1]=eigVec[i][1];
        eigVec_temp[i][2]=eigVec[i][2];
   }


   /*** make mirror images of allowed orientations ***/
   /* original eigenvectors */
   add_eigVec_List(eigVecXlist,eigVecYlist,eigVecZlist,n,eigVec_temp);


   return( 1 );
}
/* =============================================================== */


/* =============================================================== */
void  legalNotice()
{

/*
printf("\n\n |------------------------------------------------------------------|-");
printf("\n |             PALES -- Legal Notice                                  |");
printf("\n |                                                                   |");
printf("\n | This  program  and  its related software is provided  'as is'.    |");
printf("\n | It is based on a design originating in the NIH  Laboratory  of    |");
printf("\n | Chemical Physics, NIDDK.  The NIH,  the Laboratory of Chemical    |");
printf("\n | Physics,  NIDDK,  and the authors make  no warranties,  either    |");
printf("\n | express or implied,  as to any matter whatsoever with  respect    |");
printf("\n | to the software.  In particular,  any and  all  warranties  of    |");
printf("\n | merchantability and fitness for  any  particular  purpose  are    |");
printf("\n | expressly excluded.  In no event will the NIH,  the Laboratory    |");
printf("\n | of Chemical Physics,  NIDDK,  or the authors be liable for any    |");
printf("\n | loss  of  profits,  any  incidental,  special,  exemplary,  or    |");
printf("\n | consequential damages  of  any  nature  whatsoever  (including    |");
printf("\n | without  limitation,  loss  of  use  or  other  commercial  or    |");
printf("\n | research loss) arising out  of  or  relating  to  the  use  or    |");
printf("\n | performance of the software.                                      |");
printf("\n *-------------------------------------------------------------------*");
printf("\n | This program and its related software is based  on  a  design     |");
printf("\n | originating in the NIH Laboratory of Chemical Physics, NIDDK.     |");
printf("\n | It is not to be distributed or modified. (C) M. Zweckstetter 2000 |");
printf("\n *-------------------------------------------------------------------*\n\n");
printf("\n\n");
*/
FPR( stderr, "\n\n *--------------------------------------------------------------------* ");
FPR( stderr, "\n |               PALES                      (C) M. Zweckstetter 2000  | ");
FPR( stderr, "\n *--------------------------------------------------------------------*\n\n ");
FPR( stderr, "\n\n");
}
/* =============================================================== */


/* =============================================================== */
int writeEulerInfo( double irrSaupe[5], double gMag,
	double eigVal[3], double eigVec[3][3],
	double alphaList[NO_LIST], double betaList[NO_LIST],
	double gammaList[NO_LIST], char *result, int verbose )
{
   int  i,error;
   char ctemp[NAMELEN+1];

/***/
/* Initialize:
/***/

    error     = 0;
    i         = 0;
    *ctemp    = '\0';
    *result   = '\0';


/***/
/* Write eigemvalues & eigenvectors:
/***/

   SPR( result,
          "\n\nDATA IRREDUCIBLE REPRESENTATION (A0,A1R,A1I,A2R,A2I)" );

   (void) strcat( result, ctemp );


   SPR( ctemp,
          "   % .4e % .4e % .4e % .4e % .4e\n",
          irrSaupe[0], irrSaupe[1], irrSaupe[2], irrSaupe[3], irrSaupe[4] );

   (void) strcat( result, ctemp );



   SPR( ctemp,
          "DATA IRREDUCIBLE GENERAL_MAGNITUDE" );
   (void) strcat( result, ctemp );


   SPR( ctemp,
          "  % .4e\n",  gMag );
   (void) strcat( result, ctemp );


   if ( verbose )
     {

        SPR( ctemp,
          "\nREMARK Eigensystem & Euler angles for clockwise rotation about z, y', z''.\n\n" );
        (void) strcat( result, ctemp );


        SPR( ctemp,
          "DATA EIGENVALUES (Axx,Ayy,Azz)" );
        (void) strcat( result, ctemp );


        SPR( ctemp,
          "   % .4e % .4e % .4e\n",
          eigVal[2], eigVal[1], eigVal[0] );

        (void) strcat( result, ctemp );



        SPR( ctemp,
          "DATA EIGENVECTORS\n" );

        (void) strcat( result, ctemp );


        SPR( ctemp,
          "DATA EIGENVECTORS XAXIS % .4e % .4e % .4e\n",
          eigVec[0][2],
          eigVec[1][2],
          eigVec[2][2] );

        (void) strcat( result, ctemp );

        SPR( ctemp,
          "DATA EIGENVECTORS YAXIS % .4e % .4e % .4e\n",
          eigVec[0][1],
          eigVec[1][1],
          eigVec[2][1] );

        (void) strcat( result, ctemp );

        SPR( ctemp,
          "DATA EIGENVECTORS ZAXIS % .4e % .4e % .4e\n",
          eigVec[0][0],
          eigVec[1][0],
          eigVec[2][0] );

        (void) strcat( result, ctemp );


/***/
/* Write Euler angles:
/***/

        SPR( ctemp, "\n" );

        (void) strcat( result, ctemp );


        SPR( ctemp,
          "DATA Q_EULER_SOLUTIONS  %7s  %7s  %7s\n",
          "ALPHA", "BETA", "GAMMA" );

        (void) strcat( result, ctemp );

        for( i = 0; i < NO_LIST; i++ )
           {
              SPR( ctemp,
                "DATA Q_EULER_ANGLES %2d  %7.2f  %7.2f  %7.2f\n",
                i+1,
                alphaList[i],
                betaList[i],
                gammaList[i] );

              (void) strcat( result, ctemp );
           }

     } /* end of if verbose */



/***/
/* Exit point:
/***/

shutdown:

    return( error );
}
/* =============================================================== */

/* ==================   determine deviation of axes  ===================== */
int calcDcos(double eigVec [3][3], double eigVec_svd[3][3], float dcos[3])
{
  int   error;
  float x[3],y[3],z[3];
  float xb[3],yb[3],zb[3];


/***/
/* Initialization:
/***/

  error = 0;

  z[0] = eigVec[0][0];
  z[1] = eigVec[1][0];
  z[2] = eigVec[2][0];
  y[0] = eigVec[0][1];
  y[1] = eigVec[1][1];
  y[2] = eigVec[2][1];
  x[0] = eigVec[0][2];
  x[1] = eigVec[1][2];
  x[2] = eigVec[2][2];
  zb[0] = eigVec_svd[0][0];
  zb[1] = eigVec_svd[1][0];
  zb[2] = eigVec_svd[2][0];
  yb[0] = eigVec_svd[0][1];
  yb[1] = eigVec_svd[1][1];
  yb[2] = eigVec_svd[2][1];
  xb[0] = eigVec_svd[0][2];
  xb[1] = eigVec_svd[1][2];
  xb[2] = eigVec_svd[2][2];


/***/
/* Calculate directional cosines:
/***/

  dcos[0] = acos( x[0]*xb[0] + x[1]*xb[1] + x[2]*xb[2] ); 
  dcos[1] = acos( y[0]*yb[0] + y[1]*yb[1] + y[2]*yb[2] );
  dcos[2] = acos( z[0]*zb[0] + z[1]*zb[1] + z[2]*zb[2] );


/***/
/* Exit point:
/***/

shutdown:

    return( error );
}
/* =============================================================== */


/* ============   select atoms taking care of charges  ================= */
int chCoorSelect(struct PDBInfo *pdbInfo, struct PalesInfo *palesInfo, 
               int *chAtomSelCount, 
	       double *xI, double *yI, double *zI ,
               int chCount, int *chargeAtomID)
{
  int      i,id,k,error;
  int      atomCount, atomSelCount; 
  double   *xItemp, *yItemp, *zItemp;
  struct   AtomInfo *atomList;


/***/
/* Initialize
/***/

   error = 0;

   i = 0;
   k = 0;
   id = 0;
   atomCount = 0;
   atomSelCount = 0;

   atomList  = (struct AtomInfo *)NULL;

   xItemp        = (double *)NULL;
   yItemp        = (double *)NULL;
   zItemp        = (double *)NULL;


   atomCount = pdbInfo->atomCount;
   atomList  = pdbInfo->atomList;


/***/
/* Allocate temporary coordinate arrays
/***/

   xItemp       = dFltAlloc( "chCoorSelect", atomCount+1);
   yItemp       = dFltAlloc( "chCoorSelect", atomCount+1 );
   zItemp       = dFltAlloc( "chCoorSelect", atomCount+1);




   if (!xItemp || !yItemp || !zItemp )
       {
        FPR( stderr, "chCoorSelect: Error allocating coord buffers.\n" );
        error = 1;
        goto shutdown;
       }

/***/
/* Select requested coordinates (atomtype,resID,surface accessibility)
/***/

  if ( error = coorSelect( pdbInfo, palesInfo, &atomSelCount, xItemp, yItemp, zItemp) )
       {
        FPR( stderr, "Error reading selection of atoms.\n");
        goto shutdown;
       }


   if( palesInfo->elInfo.elFlag )
     {
       for( i=1; i<=chCount; i++ )
          {
		id    = chargeAtomID[i];
                xI[i] = atomList[id].x;
                yI[i] = atomList[id].y;
                zI[i] = atomList[id].z;
                chargeAtomID[i] = i;
          }

       k =  chCount;
   }
   else {
	k =  0;
   }
  
   for( i=1; i<=atomSelCount; i++)
      {
                xI[k+i] = xItemp[i];
                yI[k+i] = yItemp[i];
                zI[k+i] = zItemp[i];
      }

   (*chAtomSelCount) = atomSelCount+k;



/***/
/* Special Case
/***/

shutdown:

  (void) deAlloc( "chCoorSelect", xItemp,      sizeof(double)*(atomCount+1) );
  (void) deAlloc( "chCoorSelect", yItemp,      sizeof(double)*(atomCount+1) );
  (void) deAlloc( "chCoorSelect", zItemp,      sizeof(double)*(atomCount+1) );

    return( error );
}
/* =============================================================== */


/* ========================================================================= */
/* ========================================================================= */
/* ========================================================================= */

/* ===================== reading pdb file ============================= */
int readPDBfree( char *inName, int *atomCount, double *xI, 
	double *yI, double *zI )
/* no molecule specifiaction allowed (if more than one molecule) */
/* key word ATOM not allowed before start of coordinates */
/* Molmol:      no plus sign at amino acid type allowed */
/*              no pseudoatoms Q allowed */
{
  char  dummy[NAMELEN+1];
  int   i,error,intDummy;
  float x,y,z;
  FILE  *inUnit;


/***/
/* Initialize
/***/

  error    = 0;
  i        = 0;

  *dummy   = '\0';

  inUnit   = (FILE *) NULL;



  if ( *inName )
     {

/***/
/* Open the file with the Euler anlges alpha, beta and gamma:
/***/

       if (!(inUnit = fopen( inName, "r")))
          {
             FPR( stderr, "readPDBfree File open error: %s.\n", inName );
             error = 1;
             goto shutdown;
          }


/***/
/* Read the header:
/***/
       do
         {
            fscanf(inUnit,"%s",dummy);
         } while (strcmp(dummy,"ATOM")!=0);


/***/
/* Read atom coordinates:
/***/

       i = 0;
       while ((strcmp(dummy,"ATOM")==0)&&(i<(MAX_ATOM_NO-1)))
       { 
           i++;

           fscanf(inUnit,"%d %s %s %d %f %f %f",&intDummy, 
                   dummy,dummy,&intDummy,&x,&y,&z);


           xI[i] = x;
           yI[i] = y;
           zI[i] = z;
  
           do    /* read remaining part of line */
           {
                        fscanf(inUnit,"%s",dummy);
           }while ( (strcmp(dummy,"ATOM")!=0) && (strcmp(dummy,"END")!=0) );
       }

       if( i==(MAX_ATOM_NO-1) )
         {
            FPR( stderr, "error in readPDBfree: too many atoms.\n" );
            error = 1;
            goto shutdown;
         }
    }


  *atomCount = i;


/***/
/* Exit point:
/***/

shutdown:
 
   if (inUnit) (void) fclose( inUnit );

   return( error );
}
/* ========================================================================= */


/* ===================== read configuration file ============================ */
int   readConfig( char *inName, struct PalesInfo *palesInfo )
{
  int   error;
  float flagvalue;
  char  line[NAMELEN+1];                 /* line from input file */
  char  flagname[NAMELEN+1];
  FILE  *inUnit;


/***/
/* Initialize:
/***/

  error      = 0;
  inUnit     = (FILE *) NULL;



  if ( *inName )
     {


/***/
/* Open the file with the Euler anlges alpha, beta and gamma:
/***/

       if (!(inUnit = fopen( inName, "r")))
          {
             FPR( stderr, "readConfig File open error: %s.\n", inName );
             error = 1;
             goto shutdown;
          }


/***/
/* Read Euler anlges alpha, beta and gamma:
/***/

       while (1)
          {
                if ( fgets(line, sizeof(line), inUnit) == NULL )
                        break;

                if ( sscanf(line, "%s %f", flagname, &flagvalue) != 2 )
                   {
                      FPR( stderr, "Error: Input data in auto file have wrong format\n");
                      FPR( stderr, "Line:%s", line);
                      error = 1;
                      goto shutdown;
                   }

                if (strcmp(flagname,"CENTER_PIVOT:")== 0)          palesInfo->gridInfo.pivotCenter = flagvalue;  

 /* @@ */
else if (strcmp(flagname,"Pmax:")== 0)    Pmax = flagvalue;

        /**************************/
        /*** ensemble averaging ***/
        /* put out surface ensemble */
        else if (strcmp(flagname,"OUTPUT_ENSEMBLE:")== 0)       palesInfo->gridInfo.outEflag=flagvalue;
        /* density of dots on ensemble sphere */
        else if (strcmp(flagname,"ALIGN_ALONG_Z:")== 0)         palesInfo->gridInfo.initOrientFlag = flagvalue;
        /* density of dots on ensemble sphere */
        else if (strcmp(flagname,"FIX_POSITION_DOM1:")== 0)     palesInfo->gridInfo.fixRotDomain = flagvalue;
        /* density of dots on ensemble sphere */
        else if (strcmp(flagname,"ENSEMBLE_CLASH:")== 0)        palesInfo->gridInfo.clashFlag = flagvalue;
        /* density of dots on ensemble sphere */
        else if (strcmp(flagname,"CLASH_DISTANCE:")== 0)        palesInfo->gridInfo.clashDist = flagvalue;
        /* density of dots on ensemble sphere */
        else if (strcmp(flagname,"ENSEMBLE_DENSITY:")== 0)      palesInfo->gridInfo.densE = flagvalue;
        /* digitization of rotation about 3rd Euler angle */
        else if (strcmp(flagname,"ENSEMBLE_GAMMA:")== 0)        palesInfo->gridInfo.digPsiE=flagvalue;
        /* change rel. orientation of domains (0/1) */
        else if (strcmp(flagname,"DOMAIN_ANGLE_VAR:")== 0)      palesInfo->gridInfo.angleDomFlag=flagvalue;
        /* Euler angles for adjusting angle between 2 domains */        
        else if (strcmp(flagname,"DOM_ANG_ALPHA:")== 0)         palesInfo->gridInfo.angleDom[0] = flagvalue;
        else if (strcmp(flagname,"DOM_ANG_BETA:")== 0)          palesInfo->gridInfo.angleDom[1] = flagvalue;
        else if (strcmp(flagname,"DOM_ANG_GAMMA:")== 0)         palesInfo->gridInfo.angleDom[2] = flagvalue;
        /* origin for orientation of domains */
        else if (strcmp(flagname,"ORIGIN_DOM1:")== 0)           palesInfo->gridInfo.resDom1Origin=flagvalue;
        else if (strcmp(flagname,"ORIGIN_DOM2:")== 0)           palesInfo->gridInfo.resDom2Origin=flagvalue;
        /* domain boundaries */
        else if (strcmp(flagname,"START_RESID_1:")== 0)         palesInfo->gridInfo.startResDom1=flagvalue;
        else if (strcmp(flagname,"END_RESID_1:")== 0)           palesInfo->gridInfo.endResDom1=flagvalue;
        else if (strcmp(flagname,"START_RESID_2:")== 0)         palesInfo->gridInfo.startResDom2=flagvalue;
        else if (strcmp(flagname,"END_RESID_2:")== 0)           palesInfo->gridInfo.endResDom2=flagvalue;
        /* limits for polar angles of cone */
        else if (strcmp(flagname,"BETA_L:")== 0)        palesInfo->gridInfo.lBetaCone = flagvalue;
        else if (strcmp(flagname,"BETA_H:")== 0)        palesInfo->gridInfo.hBetaCone = flagvalue;
        else if (strcmp(flagname,"ALPHA_CENTER:")== 0)  palesInfo->gridInfo.cAlphaCone = flagvalue;
        else if (strcmp(flagname,"ALPHA_WIDTH:")== 0)   palesInfo->gridInfo.errAlphaCone = flagvalue;
        /* limits for calculating characteristic angle */
        else if (strcmp(flagname,"START_RES_ANGLE_1:")== 0)     palesInfo->gridInfo.startResPart1=flagvalue;
        else if (strcmp(flagname,"END_RES_ANGLE_1:")== 0)       palesInfo->gridInfo.endResPart1=flagvalue;
        else if (strcmp(flagname,"START_RES_ANGLE_2:")== 0)     palesInfo->gridInfo.startResPart2=flagvalue;
        else if (strcmp(flagname,"END_RES_ANGLE_2:")== 0)       palesInfo->gridInfo.endResPart2=flagvalue;

        else
           FPR( stderr, "Error: Wrong input parameter in auto file" );     
        }
     }



/***/
/* Exit point:
/***/

shutdown:
 
   if (inUnit) (void) fclose( inUnit );

   return( error );
}
/* ========================================================================= */

/* ========================================================================= */
/* ========================================================================= */
/* flexibility modules                                                       */
/* ========================================================================= */
/* ========================================================================= */


/* ============================= rotate protein =========================== */
int  rotateDomain( double angle[3], struct PDBInfo *pdbInfo,
	int resN, int resC )
{
  int     i,j;
  float   error;
  double  atom_vec[4];
  double  x_axis_tensor[4],y_axis_tensor[4],z_axis_tensor[4];
  double  Dreh_m[4][4];


  i     = 0;
  j     = 0;
  error = 0;


  /* rotate coordinate system (inverse rotation) */
  rotate_coorsystem(-1.*RAD(angle[1]),-1.*RAD(angle[2]),-1.*RAD(angle[0]),
        x_axis_tensor,y_axis_tensor,z_axis_tensor);
  for(i=0; i<3; i++)
     {
        Dreh_m[i][0] = x_axis_tensor[i+1];
        Dreh_m[i][1] = y_axis_tensor[i+1];
        Dreh_m[i][2] = z_axis_tensor[i+1];
     }


  /* calculate coordinates in rotated frame */
  for( i = 0; i < pdbInfo->atomCount; i++ )
     {
        if( pdbInfo->atomList[i].resID >= resN )
          {
             if( pdbInfo->atomList[i].resID <= resC ) 
               {
                  atom_vec[0]=pdbInfo->atomList[i].x;       /* select cartesian */
                  atom_vec[1]=pdbInfo->atomList[i].y;       /* coordinates of */
                  atom_vec[2]=pdbInfo->atomList[i].z;       /* atom i */

                  pdbInfo->atomList[i].x = 0;
                  pdbInfo->atomList[i].y = 0;
                  pdbInfo->atomList[i].z = 0;
                  for( j=0; j<3; j++ )
                     {
                        pdbInfo->atomList[i].x += Dreh_m[j][0]*atom_vec[j];
                        pdbInfo->atomList[i].y += Dreh_m[j][1]*atom_vec[j];
                        pdbInfo->atomList[i].z += Dreh_m[j][2]*atom_vec[j];
                     }
               }
             else
                 break;
          }
     }


/***/
/* Exit point:
/***/

shutdown:
 
    return( error );
}
/* ========================================================================= */


/* ============================= rotate protein =========================== */
int positionDomain(struct PDBInfo *pdbInfo, struct PalesInfo *palesInfo, int resN, int resC,
        int origin1, int origin2 )
{
  int     i,error;
  double  domain_origin[3];
  double  inertiaTensor[9];
  double  da, dr;
  double  eigVec[3][3], eigVal[3];
  static  double alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
  double  rotAngle[3];



/***/
/* Initialize:
/***/

    error       = 0;
    i           = 0;


/***/
/* Shift whole molecule to pivot point of second domain:
/***/

  for( i = 0; i < pdbInfo->atomCount; i++ )
     {
       if( pdbInfo->atomList[i].resID == origin2 )
          {
            if( strcmp(pdbInfo->atomList[i].atomName,"CA") == 0 )
              {
                 domain_origin[0] = pdbInfo->atomList[i].x;
                 domain_origin[1] = pdbInfo->atomList[i].y;
                 domain_origin[2] = pdbInfo->atomList[i].z;
                 break;
              }
          }
     }


  for( i = 0; i < pdbInfo->atomCount; i++ )
     {
        pdbInfo->atomList[i].x -= domain_origin[0];
        pdbInfo->atomList[i].y -= domain_origin[1];
        pdbInfo->atomList[i].z -= domain_origin[2];
     }



/***/
/* Calculate Inertia tensor:
/***/
 
/*
  if( error = calcInertia( pdbInfo,pdbInfo->minResID,pdbInfo->maxResID,
				inertiaTensor ) )
    {
	FPR( stderr, "Error in calcInertia: zero atoms selected." );
	goto shutdown;
    }


  if( analyzeMatrix3D(inertiaTensor, eigVec, eigVal, &da, &dr,
        alphalist, betalist, gammalist) )
    {
	FPR( stderr, "positionDomain: error in analyzeMatrix3D." );
	error = 1;
	goto shutdown;
    }
*/


 
/***/
/* Rotate selected part of molecule
/***/

  rotAngle[0]=palesInfo->gridInfo.angleDom[0];
  rotAngle[1]=palesInfo->gridInfo.angleDom[1];
  rotAngle[2]=palesInfo->gridInfo.angleDom[2];
  if( error = rotateDomain(rotAngle,pdbInfo,resN,resC) )
    {
	FPR( stderr, "initOrient: error in analyzeMatrix3D." );
	goto shutdown;
    }



/***/
/* Put protein in original translational postion:
/***/

  for( i = 0; i < pdbInfo->atomCount; i++ )
     {
       if( pdbInfo->atomList[i].resID == origin1 )
          {
            if( strcmp(pdbInfo->atomList[i].atomName,"CA") == 0 )
              {
                 domain_origin[0] = pdbInfo->atomList[i].x;
                 domain_origin[1] = pdbInfo->atomList[i].y;
                 domain_origin[2] = pdbInfo->atomList[i].z;
                 break;
              }
          }
     }


  for( i = 0; i < pdbInfo->atomCount; i++ )
     {
        pdbInfo->atomList[i].x -= domain_origin[0];
        pdbInfo->atomList[i].y -= domain_origin[1];
        pdbInfo->atomList[i].z -= domain_origin[2];
     }


/***/
/* Exit point:
/***/

shutdown:
 
    return( error );
}
/* ========================================================================= */

/* ========================================================================= */
int simulateAlignment( struct PalesInfo *palesInfo, int atomSelCount, 
        double *xI, double *yI, double *zI , 
	int chCount, int *chargeAtomID, float *charges,
	int dotCount, double *ddots,
	int digPotential, double *potential, double offset, double stepsize,
	double saupe[5] )
{
   int    i, error, density;
   float  center_of_mass[3], maxRadius, minRadius, distance;
   double radius,gridSpacing,wvAmount,molarRatio;


/***/
/* Initialize:
/***/

    error   = 0;
    i       = 0;
    density = 0;
    distance = 0;
    maxRadius = 0.0;
    minRadius = 100000.0;
    for( i=0; i<3; i++ )     center_of_mass[i]=0;
 

/***/
/* Shift protein to geometric center
/* Determine maximum and minimum radius of protein
/***/

   /*** calculate geometric center of positioning helix ***/
   for( i=0; i<3; i++ )     center_of_mass[i]=0;
   for( i=1; i<=atomSelCount; i++ )
      {
                center_of_mass[0] += xI[i];
                center_of_mass[1] += yI[i];
                center_of_mass[2] += zI[i];
      }

   center_of_mass[0] /= atomSelCount;
   center_of_mass[1] /= atomSelCount;
   center_of_mass[2] /= atomSelCount;


  /*** shift selected part of protein to geometric center ***/
  for( i=1; i<=atomSelCount; i++)
     {
                xI[i] -= center_of_mass[0];
                yI[i] -= center_of_mass[1];
                zI[i] -= center_of_mass[2];
     }


  /*** determine extrema of length of protein ***/
  minRadius = 100000.0;
  maxRadius = 0.0;
  for( i=1; i<=atomSelCount; i++ )
     {
        distance = xI[i]*xI[i]+yI[i]*yI[i]+zI[i]*zI[i];
        if(distance > maxRadius) maxRadius = distance;
        if(distance < minRadius) minRadius = distance;
     }
  minRadius = sqrt(minRadius);
  maxRadius = sqrt(maxRadius);



/***/
/* Sweep center of mass through 1D lattice
/***/
 
         radius      = palesInfo->radiusAtom;
         gridSpacing = palesInfo->gridSpacing;
         wvAmount    = palesInfo->concentration;
         molarRatio  = palesInfo->molarRatio;


         if( error = stericAlignment( palesInfo, palesInfo->orderDiscRim,palesInfo->lcType,palesInfo->rSingleCell,palesInfo->radiusModel,
                    radius,gridSpacing,wvAmount,molarRatio,
	            atomSelCount,xI,yI,zI,maxRadius+palesInfo->radiusAtom,minRadius+palesInfo->radiusAtom,
                    palesInfo->digPsiNSC,dotCount,ddots,saupe ) )
           {
	     FPR( stderr, "simulateAlignment: error in stericAlignment." );
	     goto shutdown;
           } 

/*
         ssia( palesInfo,atomSelCount,xI,yI,zI,
               maxRadius+palesInfo->radiusAtom,minRadius+palesInfo->radiusAtom,
               dotCount,ddots,saupe );
*/



    /* scale Saupe matrix according to bicelle order */
    for( i=0; i<5; i++ )    saupe[i] *= palesInfo->lcOrder; 


/***/
/* Exit point:
/***/

shutdown:
 
    return( error );
}
/* ========================================================================= */


/* ========================================================================= */
int writeTensor( outDName, access, dTab, remark, verbose )

   struct DTable  *dTab;
   char   *outDName, *remark, *access;
   int    verbose;
{
    struct DInfo *dObsList, *dCalcList, *dObsPtr, *dCalcPtr;
    float  shift, d, rSumD, pSumD;
    int    i, nD;

    FILE   *outUnit;

/***/
/* Calculate statistics of fit.
/***/

    if (!dTab)            return( 0 );
    if (!dTab->dCount)    return( 0 );

    if (!dTab->dObsList)  return( 0 );
    if (!dTab->dCalcList) return( 0 );

    dObsList  = dTab->dObsList;
    dCalcList = dTab->dCalcList;

    (void) getStatDC( dTab->dObsList,
                      dTab->dCalcList,
                      dTab->dCount,
                      &rSumD, &pSumD, &nD );

/***/
/* Report of simulated couplings:
/***/

    outUnit = stdout;

    if (*outDName)
       {
        if (!(outUnit = fopen( outDName, access )))
           {
            FPR( stderr, "DC Error opening output file %s\n", outDName );
            outUnit = stdout;
           }
       }
 
#ifdef DC_NEW_CONST
    FPR( outUnit, "REMARK Standard Dipolar Coupling Simulation.\n" );
#else
    FPR( outUnit, "REMARK Dipolar Coupling Simulation.\n" );
#endif

    FPR( outUnit, "REMARK %s\n\n", remark );

    FPR( outUnit, "DATA N    %d\n",       nD );
    FPR( outUnit, "DATA RMS  %.3f\n",     rSumD );
    FPR( outUnit, "DATA Chi2 %.3f\n\n",   pSumD );

    switch( dTab->dTensor.tensorMode )
       {
        case DC_STERIC:
           FPR( outUnit, "DATA TENSOR_MODE STERIC" );
           FPR( outUnit, " (Zweckstetter Method)\n" );
           break;

        case DC_ELECTROSTATIC:
           FPR( outUnit, "DATA TENSOR_MODE ELECTROSTATIC" );
           FPR( outUnit, " (Zweckstetter Method)\n" );
           break;

        case DC_SVD:
           FPR( outUnit, "DATA TENSOR_MODE SVD (Order Matrix Method)\n" );
           break;

        case DC_SAUPE:
           FPR( outUnit, "DATA TENSOR_MODE SAUPE" );
           FPR( outUnit, " (Fixed Saupe Values)\n" );
           break;

        case DC_FIXED_ALL:
           FPR( outUnit, "DATA TENSOR_MODE FIXED" );
           FPR( outUnit, " (Fixed Da and Dr, Fixed Rotations)\n" );
           break;

        case DC_FIXED_DADR:
           FPR( outUnit, "DATA TENSOR_MODE DADR" );
           FPR( outUnit, " (Fixed Da and Dr, Variable Rotations)\n\n" );

           FPR( outUnit,
                "DATA PSI_RANGE   % 8.2f % 8.2f\n",
                dTab->dTensor.psiMin, dTab->dTensor.psiMax );

           FPR( outUnit,
                "DATA THETA_RANGE % 8.2f % 8.2f\n",
                dTab->dTensor.thetaMin, dTab->dTensor.thetaMax );

           FPR( outUnit,
                "DATA PHI_RANGE   % 8.2f % 8.2f\n",
                dTab->dTensor.phiMin, dTab->dTensor.phiMax );
           break;

        case DC_FIXED_NONE:
           FPR( outUnit, "DATA TENSOR_MODE NOFIXED" );
           FPR( outUnit, " (Variable Da and Dr, Variable Rotations)\n\n" );

           FPR( outUnit,
                "DATA DA_RANGE    % .4e % .4e\n",
                dTab->dTensor.daMin, dTab->dTensor.daMax );

           FPR( outUnit,
                "DATA DR_RANGE    % .4e % .4e\n",
                dTab->dTensor.drMin, dTab->dTensor.drMax );

           FPR( outUnit,
                "DATA PSI_RANGE   % 8.2f % 8.2f\n",
                dTab->dTensor.psiMin, dTab->dTensor.psiMax );

           FPR( outUnit,
                "DATA THETA_RANGE % 8.2f % 8.2f\n",
                dTab->dTensor.thetaMin, dTab->dTensor.thetaMax );

           FPR( outUnit,
                "DATA PHI_RANGE   % 8.2f % 8.2f\n",
                dTab->dTensor.phiMin, dTab->dTensor.phiMax );
           break;

        default:
           break;
       }

    FPR( outUnit, "\n" );

    FPR( outUnit,
         "DATA SAUPE  % .4e % .4e % .4e % .4e % .4e\n\n",
         dTab->dTensor.saupe[0], 
         dTab->dTensor.saupe[1], 
         dTab->dTensor.saupe[2], 
         dTab->dTensor.saupe[3], 
         dTab->dTensor.saupe[4] );


    (void) reportEuler( outUnit, &dTab->dTensor );

    FPR( outUnit, "VARS   " );

    PrintAtomIVars( outUnit, dTab );
    PrintAtomJVars( outUnit, dTab );

    if (verbose)
       FPR( stderr, " DI D_OBS D D_DIFF DD W\n" );
    else
       FPR( stderr, " DI D DD W\n" );

    FPR( outUnit, "FORMAT " );
   
    PrintAtomIFmt( outUnit, dTab );
    PrintAtomJFmt( outUnit, dTab );

    if (verbose)
       FPR( stderr, "%%9.2f %%9.3f %%9.3f %%9.3f %%.2f %%.2f\n" );
    else
       FPR( stderr, "%%9.2f %%9.3f %%9.3f %%.2f\n" );

    dObsPtr  = dObsList;
    dCalcPtr = dCalcList;
 
    for( i = 0; i < dTab->dCount; i++ )
       {
        if (dCalcPtr->valid)
           {
            if (dObsPtr->d < 999.0 && dCalcPtr->d < 999.0)
               d = dObsPtr->d - dCalcPtr->d;
            else
               d = DC_UVAL;

            PrintAtomIVals( outUnit, dObsPtr, dTab );
            PrintAtomJVals( outUnit, dObsPtr, dTab );

            if (verbose)
               { 
                FPR( outUnit,
                     "% 9.2f % 9.4f % 9.4f % 9.4f %7.4f %.2f\n",
                     dCalcPtr->di,
                     dObsPtr->d,
                     dCalcPtr->d,
                     d,
                     dObsPtr->sigma,
                     dObsPtr->w );
               }
            else
               {
                FPR( outUnit,
                     " % 9.2f % 9.4f %7.4f %.2f\n",
                     dCalcPtr->di,
                     dCalcPtr->d,
                     dObsPtr->sigma,
                     dObsPtr->w );
               }
           }

        dObsPtr++; dCalcPtr++;
       }

    if (outUnit != stdout) (void) fclose( outUnit );

    return( 0 );
}
/* ========================================================================= */


/* ========================================================================= */
int getparams( argc, argv,
               inpName, potName, pkaName, inEName, flexName, oDcEName,
               pdbFName, pRotName, outEName, outEAxName, aveName,
               inDName, inCSName, outDName, outCSName, outSName,
               surfPath, atomNames, resNames, pdbList, pdbPath, pdbName,
               refName, dTensor, fitFlag, wList, wCount, csThresh, validThresh,
               mapCount, outMapName, userSaupe1, userSaupe2,
               outAName, remark, atomNameS1, atomNameS2,
               inD1Name, inD2Name, palesInfo, sInfo, dcInfo )

   int   argc, *fitFlag, *wCount, *mapCount;

   char  *inpName, *potName, *pkaName, *inEName, *flexName, *oDcEName,
         *pdbFName, *pRotName, *outEAxName, *aveName, *outEName;

   char  **remark;

   char  **argv, *inDName, *inCSName, *outDName, *outCSName, *outSName,
         *surfPath, *atomNames, *resNames, *pdbList, *pdbPath, *pdbName,
         *refName, *outMapName, *outAName, *atomNameS1, *atomNameS2,
         *inD1Name, *inD2Name;

   double  userSaupe1[5], userSaupe2[5];

   float **wList, *csThresh, *validThresh;
   
   struct DTensor *dTensor;
   struct PDBSelectInfo *sInfo;
   struct PalesInfo *palesInfo;
   struct DcInfo     *dcInfo;
{
    int i, loc1, loc2, loc, nextLoc, selectCount;
    char *sPtr;

float dummy;

/***/
/* Install defaults;
/***/

    selectCount          = 0;
    *mapCount = 0;




    if (argc == 1)
       {
         (void) printInfo();

         return( 1 );
       }

    if (flagLoc( argc, argv, "-help" ))
       {
         (void) printInfo();

         return( 1 );
       }


    palesInfo->palesMode = PALES_STERIC;

    if (flagLoc( argc, argv, "-stPales" ))
       {
          palesInfo->palesMode = PALES_STERIC;
       }

    if (flagLoc( argc, argv, "-stPalesFree" ))
       {
          palesInfo->palesMode = PALES_STERIC_FREE;
       }

    if (flagLoc( argc, argv, "-bestFit" ))
       {
          palesInfo->palesMode = PALES_DC;
       }
    
    if (flagLoc( argc, argv, "-anA" ))
       {
          palesInfo->palesMode = PALES_SAUPE_CP;
          palesInfo->c1        = 1.0;
          palesInfo->c2        = 1.0;
          palesInfo->opCode    = OP_SUB;
       }

    if (flagLoc( argc, argv, "-anDC" ))
       {
          palesInfo->palesMode = PALES_DC_CP;
          palesInfo->c1        = 0.0;
          palesInfo->c2        = 1.0;
          palesInfo->opCode    = OP_ADD;
       }


    (void) intArgD( argc, argv, "-map",    mapCount );

    if (flagLoc( argc, argv, "-pdb" ))
       {
        dTensor->eulerFlag = 1;
        *fitFlag           = 1;
       }
    else
       {
        dTensor->eulerFlag = 0;
        *fitFlag           = 0;
       }


    (void) getPDBSelectArgs( argc, argv, sInfo );

 
    dTensor->eulerFlag = 1;

    (void) getDCTensorParams( argc, argv, dTensor );

    palesInfo->eulerFlag = dTensor->eulerFlag;


    if (loc = flagLoc( argc, argv, "-dbg" ))
       {
        if (loc + 1 == nextFlag( argc, argv, loc ))
           dbgFlag = 1;
        else
           (void) intArgD( argc, argv, "-dbg", &dbgFlag );
       }

    if (loc = flagLoc( argc, argv, "-noverb" ))
       {
        if (loc + 1 == nextFlag( argc, argv, loc ))
           palesInfo->verbFlag = 0;
        else
           (void) intArgD( argc, argv, "-noverb", &palesInfo->verbFlag );
       }

    if (sPtr = getenv( "PDBH_DIR" ))
       {
        (void) strcpy( pdbPath, sPtr );
       }

    if (sPtr = getenv( "PDBH_LIST" ))
       {
        (void) strcpy( pdbList, sPtr );
       }

    if (sPtr = getenv( "SURF_DIR" ))
       {
        (void) strcpy( surfPath, sPtr );
       }

    *wList = rListArg( argc, argv, "-w", wCount );

    (void) strArgD( argc, argv, "-in",     inDName );
    (void) strArgD( argc, argv, "-inD",    inDName );
    (void) strArgD( argc, argv, "-inCS",   inCSName );
    (void) strArgD( argc, argv, "-out",    outDName );
    (void) strArgD( argc, argv, "-outD",   outDName );
    (void) strArgD( argc, argv, "-outS",   outSName );
    (void) strArgD( argc, argv, "-outCS",  outCSName );
    (void) strArgD( argc, argv, "-pdb",    pdbName );
    (void) strArgD( argc, argv, "-ref",    refName );
    (void) strArgD( argc, argv, "-path",   pdbPath );
    (void) strArgD( argc, argv, "-list",   pdbList );
    (void) strArgD( argc, argv, "-sPath",  surfPath );
    (void) strArgD( argc, argv, "-sList",  atomNames );
    (void) strArgD( argc, argv, "-rList",  resNames );

    if (!strcasecmp( atomNames, "All" ))
       {
        (void) strcpy( atomNames, "HA CA CB C N" );
       }

    if (!strcasecmp( resNames, "All" ))
       {
        (void) strcpy( resNames, "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE " );
        (void) strcat( resNames, "LEU LYS MET PHE PRO SER THR TRP TYR VAL " );
        (void) strcat( resNames, "xPRO" );
       }

    (void) intArgD( argc, argv, "-dStep", &dcInfo->segStep );
    (void) intArgD( argc, argv, "-sStep", &dcInfo->segStep );
    (void) intArgD( argc, argv, "-d1",    &dcInfo->segFirst );
    (void) intArgD( argc, argv, "-dN",    &dcInfo->segLast );
    (void) intArgD( argc, argv, "-s1",    &dcInfo->segFirst );
    (void) intArgD( argc, argv, "-sN",    &dcInfo->segLast );
    (void) intArgD( argc, argv, "-sID",   &dcInfo->segID );
    (void) intArgD( argc, argv, "-a1",    &dcInfo->aFirst );

    (void) fltArgD( argc, argv, "-csThresh", csThresh );
    (void) fltArgD( argc, argv, "-vcThresh", validThresh );


    if (flagLoc( argc, argv, "-excl" ))      dcInfo->excludeMode = DC_EXCLUDE;
    if (flagLoc( argc, argv, "-inc" ))       dcInfo->excludeMode = DC_INCLUDE_ALL;

    if (flagLoc( argc, argv, "-fit" ))       *fitFlag   = 1;
    if (flagLoc( argc, argv, "-nofit" ))     *fitFlag   = 0;

    (void) ptrArgD( argc, argv, "-rem", remark );
 

/***/
/* PALES options
/***/

    (void) strArgD( argc, argv, "-inp",    inpName );
    (void) strArgD( argc, argv, "-outD",   outDName );
    (void) strArgD( argc, argv, "-outMap", outMapName );
    (void) strArgD( argc, argv, "-pdbF",   pdbFName);
    (void) strArgD( argc, argv, "-pdbRot", pRotName);

    (void) strArgD( argc, argv, "-aD1",    atomNameS1 );
    (void) strArgD( argc, argv, "-aD2",    atomNameS2 );


    (void) strArgD( argc, argv, "-inD1",   inD1Name );
    (void) strArgD( argc, argv, "-inD2",   inD2Name );
    (void) strArgD( argc, argv, "-outD",   outDName );


    if (flagLoc( argc, argv, "-bic" ))
      {
        palesInfo->lcType      = 1;
        palesInfo->radiusModel = 20.0;
      }
    if (flagLoc( argc, argv, "-pf1" ))
      {
        palesInfo->lcType      = 2;
        palesInfo->radiusModel = 33.5;
      }

    (void) fltArgD( argc, argv, "-wv",       &palesInfo->concentration );
    (void) intArgD( argc, argv, "-dot",      &palesInfo->densityNSC );
    (void) intArgD( argc, argv, "-digPsi",   &palesInfo->digPsiNSC );
    (void) fltArgD( argc, argv, "-dGrid",    &palesInfo->gridSpacing );
    (void) fltArgD( argc, argv, "-rM",       &palesInfo->radiusModel );
    (void) fltArgD( argc, argv, "-lcS",      &palesInfo->lcOrder );

    if (flagLoc( argc, argv, "-H" ))         palesInfo->Hflag    = 1;
    if (flagLoc( argc, argv, "-noH" ))       palesInfo->Hflag    = 0;
    if (flagLoc( argc, argv, "-surf" ))      palesInfo->surfFlag = 1;
    if (flagLoc( argc, argv, "-nosurf" ))    palesInfo->surfFlag = 0;

    (void) fltArgD( argc, argv, "-rA",       &palesInfo->radiusAtom );

    (void) fltArgD( argc, argv, "-rim",      &palesInfo->includeRim );
    (void) fltArgD( argc, argv, "-pRatio",   &palesInfo->molarRatio );
    (void) fltArgD( argc, argv, "-lM",       &palesInfo->lengthModel );
    (void) fltArgD( argc, argv, "-lM",       &palesInfo->lengthModel );
    (void) fltArgD( argc, argv, "-rimDisc",  &palesInfo->orderDiscRim );



/***/
/* Extract saupe values for camparing order matrices:
/***/

    loc  = 0;
    loc1 = flagLoc( argc, argv, "-inS1" );

    if (loc1) loc = loc1;

    if (loc) 
       {
        nextLoc = nextFlag( argc, argv, loc );

        if (nextLoc - loc != 6)
           {
            FPR( stderr, "DYNAMO Error: -inS1 option requires 5 values.\n" );
            return( 1 );
           }

        for( i = 0; i < 5; i++ )
           {
            userSaupe1[i] = ATOF( getNthArg( argc, argv, loc+i+1 ));
           }
       }



    loc  = 0;
    loc2 = flagLoc( argc, argv, "-inS2" );

    if (loc2) loc = loc2;

    if (loc) 
       {
        nextLoc = nextFlag( argc, argv, loc );

        if (nextLoc - loc != 6)
           {
            FPR( stderr, "DYNAMO Error: -inS2 option requires 5 values.\n" );
            return( 1 );
           }

        for( i = 0; i < 5; i++ )
           {
            userSaupe2[i] = ATOF( getNthArg( argc, argv, loc+i+1 ));
           }
       }

    (void) strArgD( argc, argv, "-outA", outAName );


    if (flagLoc( argc, argv, "-add" )) palesInfo->opCode = OP_ADD;
    if (flagLoc( argc, argv, "-sub" )) palesInfo->opCode = OP_SUB;
    if (flagLoc( argc, argv, "-mul" )) palesInfo->opCode = OP_MUL;
    if (flagLoc( argc, argv, "-div" )) palesInfo->opCode = OP_DIV;
    if (flagLoc( argc, argv, "-max" )) palesInfo->opCode = OP_MAX;

    (void) fltArgD( argc, argv, "-c1", &palesInfo->c1 );
    (void) fltArgD( argc, argv, "-c2", &palesInfo->c2 );


    return( 0 );
}
/* ========================================================================= */

/* ========================================================================= */
int printInfo( )
{

    FPR( stderr, "\nPALES Prediction of ALignmEnt from Structure.\n\n" );

    FPR( stderr, "\nARGUMENTS FOR STERIC PALES\n" );
    FPR( stderr, " -stPales          [on]      Steric PALES Mode.\n" );
    FPR( stderr, " -pdb     pName    [None]    PDB Input (Required).\n" );
    FPR( stderr, " -inD     dName    [None]    Dipolar Coupling (Input).\n" );
    FPR( stderr, " -inp     inpName  [None]    Configuration file (Input).\n" );
    FPR( stderr, " -outD    dName    [None]    Dipolar Coupling (Output).\n" );
    FPR( stderr, " -pdbRot  pRotName [None]    PDB rotated (Output).\n" );
    FPR( stderr, "Simulation Parameters:\n" );
    FPR( stderr, " -bic              [on]      Wall LC Model.\n" );
    FPR( stderr, " -pf1              [off]     Rod LC Model.\n" );
    FPR( stderr, " -wv      conc     [0.05]    Liquid Crystal Concentration [mg/ml].\n" );
    FPR( stderr, " -dot     density  [100]     Orientation Density (Sphere).\n" );
    FPR( stderr, " -digPsi  digPsi   [18]      Orientation Density (3rd angle).\n" );
    FPR( stderr, " -dGrid   gridSpac [0.2]     Grid Spacing [A].\n" );
    FPR( stderr, " -rM      radius   [20.0]    Half Model Thickness [A].\n" );
    FPR( stderr, " -lcS     order    [0.8]     Liquid Crystal Order [0-1].\n" );
    FPR( stderr, " -surf             [on]      Select Surface Accessible Atoms.\n" );
    FPR( stderr, " -nosurf           [off]     Don't Select Surface Accessible Atoms.\n" );
    FPR( stderr, " -rA      radius   [0.0]     Atom radius [A].\n" );
    FPR( stderr, "PDB Selection Flags:\n" );
    FPR( stderr, " -H                [off]     Simple PDB Selection (all Atoms).\n" );
    FPR( stderr, " -noH              [on]      Simple PDB Selection (all Atoms except H).\n" );
    FPR( stderr, " -res     resID              Single Residue for Selection.\n" );
    FPR( stderr, " -r1      res1     [0]       First Residue in Segment.\n" );
    FPR( stderr, " -rN      resN     [0]       Last Residue in Segment.\n" );
    FPR( stderr, " -atom    aList    [*]       Atom Name List.\n" );
    FPR( stderr, " -resName rList    [*]       Res Name List.\n" );
    FPR( stderr, " -chain   cList    [*]       Chain Name List.\n" );
    FPR( stderr, " -seg     sList    [*]       Seg Name List.\n" );
    FPR( stderr, " -id      sID      [*]       Selection ID.\n" );
    FPR( stderr, " -inv                        Invert Selection.\n" );
    FPR( stderr, "PDB Special Selection; No Other Selection Performed:\n" );
    FPR( stderr, " -clear                      Clear Selection for All Atoms.\n" );
    FPR( stderr, " -set                        Set Selection for All Atoms.\n" );
    FPR( stderr, "PDB Name Matching Modes:\n" );
    FPR( stderr, " -exact                      Match Names Exactly (Default).\n" );
    FPR( stderr, " -regex                      Regular Expression Pattern Matching.\n" );
    FPR( stderr, " -glob                       Glob-style Pattern Matching.\n" );
    FPR( stderr, " -aExact                     Match Atom Names Exactly.\n" );
    FPR( stderr, " -aRegex                     Atom Regular Expression Pattern Matching.\n" );
    FPR( stderr, " -aGlob                      Atom Glob-style Pattern Matching.\n" );
    FPR( stderr, " -rExact                     Match Res Names Exactly.\n" );
    FPR( stderr, " -rRegex                     Res Regular Expression Pattern Matching.\n" );
    FPR( stderr, " -rGlob                      Res Glob-style Pattern Matching.\n" );
    FPR( stderr, " -cExact                     Match Chain Names Exactly.\n" );
    FPR( stderr, " -cRegex                     Chain Regular Expression Pattern Matching.\n" );
    FPR( stderr, " -cGlob                      Chain Glob-style Pattern Matching.\n" );
    FPR( stderr, " -sExact                     Match Seg Names Exactly.\n" );
    FPR( stderr, " -sRegex                     Seg Regular Expression Pattern Matching.\n" );
    FPR( stderr, " -sGlob                      Seg Glob-style Pattern Matching.\n" );
    FPR( stderr, "DC Selection Flags:\n" );
    FPR( stderr, " -a1      aRes1    [0]       First Res in PDB.\n" );
    FPR( stderr, " -s1      res1     [0]       First Res in Segment.\n" );
    FPR( stderr, " -sN      resN     [0]       Last Res in Segment.\n" );
    FPR( stderr, " -aD1     atom     {All}     DC Atom Type 1 used.\n" );
    FPR( stderr, " -aD2     atom     {All}     DC Atom Type 2 used.\n" );


    FPR( stderr, "\n\nARGUMENTS FOR FREE FORMAT, STERIC PALES\n" );
    FPR( stderr, " -stPalesFree      [off]     Free Format Steric PALES Mode.\n" );
    FPR( stderr, " -pdbF    pdbFName [None]    PDB Input (free format) (Required).\n" );
    FPR( stderr, " -outA    outAName [None]    Predicted Saupe Matrix (Output).\n" );
    FPR( stderr, "Simulation Parameters:\n" );
    FPR( stderr, " -bic              [on]      Wall LC Model.\n" );
    FPR( stderr, " -pf1              [off]     Rod LC Model.\n" );
    FPR( stderr, " -wv      conc     [0.05]    Liquid Crystal Concentration [mg/ml].\n" );
    FPR( stderr, " -dot     density  [100]     Orientation Density (Sphere).\n" );
    FPR( stderr, " -digPsi  digPsi   [18]      Orientation Density (3rd angle).\n" );
    FPR( stderr, " -dGrid   gridSpac [0.2]     Grid Spacing [A].\n" );
    FPR( stderr, " -rM      radius   [20.0]    Half Model Thickness [A].\n" );
    FPR( stderr, " -lcS     order    [0.8]     Liquid Crystal Order [0-1].\n" );
    FPR( stderr, " -surf             [on]      Select Surface Accessible Atoms.\n" );
    FPR( stderr, " -nosurf           [off]     Don't Select Surface Accessible Atoms.\n" );
    FPR( stderr, " -rA      radius   [0.0]     Atom radius [A].\n" );


    FPR( stderr, "\n\nARGUMENTS FOR BEST-FIT (SVD/POWELL)\n" );
    FPR( stderr, " -bestFit        [off]       Pure Best-Fit Mode.\n" );
    FPR( stderr, " -pdb     pName  [None]      PDB Input (Required).\n" );
    FPR( stderr, " -inD     dName  [None]      Dipolar Coupling (Input).\n" );
    FPR( stderr, " -outD    dName  [None]      Best-fit Coupling (Output).\n" );
    FPR( stderr, " -map     mapCount  [0]      Mapping of Error Space.\n" );
    FPR( stderr, " -outMap  mName  [None]      World Map Coordinates (Output).\n" );
    FPR( stderr, "Alignment Tensor Specification Methods:\n" );
    FPR( stderr, " -svd         SVD Calc of Saupe Matrix (Default).\n" );
    FPR( stderr, " -saupe sList List of 5 Saupe Values (also: -align).\n" );
    FPR( stderr, " -fixed       Fixed Da, Dr, and Rotations.\n" );
    FPR( stderr, " -nofixed     Variable Da, Dr, and Rotations.\n" );
    FPR( stderr, " -dadr        Fixed Da and Dr, Variable Rotations.\n" );
    FPR( stderr, "Alignment Tensor Magnitude and Rhombicity Limits:\n" );
    FPR( stderr, " -da          [1.0e-3]      Initial DA Value.\n" );
    FPR( stderr, " -daMin      [-1.0e-2]      Minimum Initial DA Value.\n" );
    FPR( stderr, " -daMax       [1.0e-2]      Maximum Initial DA Value.\n" );
    FPR( stderr, " -dr          [1.0e-4]      Initial DR Value.\n" );
    FPR( stderr, " -drMin      [-1.0e-2]      Minimum Initial DR Value.\n" );
    FPR( stderr, " -drMax       [1.0e-2]      Maximum Initial DR Value.\n" );
    FPR( stderr, "Rotation Angle Starting Values and Limits:\n" );
    FPR( stderr, " -psi           [45.0]      Initial X-Axis Rotation.\n" );
    FPR( stderr, " -psiMax       [360.0]      Maximum X-Axis Rotation.\n" );
    FPR( stderr, " -psiMin      [-360.0]      Minimum X-Axis Rotation.\n" );
    FPR( stderr, " -theta         [45.0]      Initial Y-Axis Rotation.\n" );
    FPR( stderr, " -thetaMax     [360.0]      Maximum Y-Axis Rotation.\n" );
    FPR( stderr, " -thetaMin    [-360.0]      Minimum Y-Axis Rotation.\n" );
    FPR( stderr, " -phi           [45.0]      Initial Z-Axis Rotation.\n" );
    FPR( stderr, " -phiMax       [360.0]      Maximum Z-Axis Rotation.\n" );
    FPR( stderr, " -phiMin      [-360.0]      Minimum Z-Axis Rotation.\n" );
    FPR( stderr, "Selection Flags:\n" );
    FPR( stderr, " -a1      aRes1  [0]         First Res in PDB.\n" );
    FPR( stderr, " -s1      res1   [0]         First Res in Segment.\n" );
    FPR( stderr, " -sN      resN   [0]         Last Res in Segment.\n" );
    FPR( stderr, " -aD1     atom   {All}       DC Atom Type 1 used.\n" );
    FPR( stderr, " -aD2     atom   {All}       DC Atom Type 2 used.\n" );


    FPR( stderr, "\n\nARGUMENTS FOR ANALYSIS OF ORDER MATRIX\n" );
    FPR( stderr, " -anA        [off]           Analysis of 2 Saupe matrices.\n" );
    FPR( stderr, " -inS1 sList List of 5 Saupe Values (Input).\n" );
    FPR( stderr, " -inS2 sList List of 5 Saupe Values (Input).\n" );
    FPR( stderr, " -outA       outAName [None] Saupe Analysis (Output).\n" );
    FPR( stderr, "Action Flags:\n" );
    FPR( stderr, " -add   Action: inS1 + inS2\n" );
    FPR( stderr, " -sub   Action: inS1 - inS2 (Default)\n" );
    FPR( stderr, " -mul   Action: inS1 * inS2\n" );
    FPR( stderr, " -div   Action: inS1 / inS2\n" );
    FPR( stderr, " -max   Action: Max(inS1, inS2)\n" );
    FPR( stderr, "            (Tests general magnitude)\n" );
    FPR( stderr, "Other Options:\n" );
    FPR( stderr, " -c1    scale [1.0]  Scale factor for inS1\n" );
    FPR( stderr, " -c2    scale [1.0]  Scale factor for inS2\n" );


    FPR( stderr, "\n\nARGUMENTS FOR ANALYSIS OF DIPOLAR COUPLINGS\n" );
    FPR( stderr, " -anDC       [off]           Analysis of 2 Sets of Dipolar Couplings.\n" );
    FPR( stderr, " -inD1       d1Name   [None] Dipolar Coupling Set 1 (Input).\n" );
    FPR( stderr, " -inD2       d2Name   [None] Dipolar Coupling Set 2 (Input).\n" );
    FPR( stderr, " -outD       outDName [None] DC Analysis (Output).\n" );
    FPR( stderr, "Selection Flags:\n" );
    FPR( stderr, " -s1      res1     [0]       First Res in Segment.\n" );
    FPR( stderr, " -sN      resN     [0]       Last Res in Segment.\n" );
    FPR( stderr, " -aD1     atom     {All}     DC Atom Type 1 used.\n" );
    FPR( stderr, " -aD2     atom     {All}     DC Atom Type 2 used.\n" );
    FPR( stderr, "Action Flags:\n" );
    FPR( stderr, " -add   Action: d1Name + d2Name (Default)\n" );
    FPR( stderr, " -sub   Action: d1Name - d2Name\n" );
    FPR( stderr, " -mul   Action: d1Name * d2Name\n" );
    FPR( stderr, " -div   Action: d1Name / d2Name\n" );
    FPR( stderr, "Other Options:\n" );
    FPR( stderr, " -c1    scale [0.0]  Scale factor for d1Name\n" );
    FPR( stderr, " -c2    scale [1.0]  Scale factor for d2Name\n" );


    FPR( stderr, "\n\nCOMMOM ARGUMENTS\n" );
    FPR( stderr, " -euler    Solve for Euler Angles, Da, Dr (Default).\n" );
    FPR( stderr, " -noeuler  Don't Solve for Euler Angles.\n" );
    FPR( stderr, " -w     wList              List of Weights for Segment.\n" );
    FPR( stderr, " -excl                     Exclude Outside Range.\n" );
    FPR( stderr, " -incl                     Include Any Within Range.\n" );
    FPR( stderr, " -noverb                   Verbose mode OFF.\n" );
    FPR( stderr, " -rem   remTxt             PDB Remark Text.\n" );

    return( 0 );
}
/* ========================================================================= */


/* ========================================================================= */
int simDC( struct PDBInfo *pdbInfo, struct DTable *dTab,
           struct PalesInfo *palesInfo, char *remark,
           char *outDName, char *pdbRotName, int verbose )
{
   char    ctemp[6000], result[6000];
   int     i,n,error;
   float   rSum, pSum;
   double  rotAngle[3];
   double  eigVal[3], eigVec[3][3];
   double  gMag, da, dr, saupe[5], irrSaupe[5];
   static  double alphaList[NO_LIST], betaList[NO_LIST], gammaList[NO_LIST];
 

/***/
/* Initialize:
/***/

   i         = 0;
   n         = 0;
   error     = 0;
   rSum      = 0;
   pSum      = 0;
   da        = 0;
   dr        = 0;

   *ctemp    = '\0';
   *result   = '\0';



/***/
/* Best-fit dipolar couplings:
/***/

   (void) fitD( pdbInfo, dTab, verbose );



/***/
/* Analyze predicted dipolar couplings:
/***/

/***/
/* Write results:
/***/

    (void) writeDCmod( outDName, "w", dTab, palesInfo, 
                    result );



/***/
/* Write molecule rotated such that alignment tensor coincide with laboratory frame:
/***/

   if( error = storePDBCoords( pdbInfo ) )
     {
        FPR( stderr, "matrix: error in storePDBCoords.\n" );
        goto shutdown;
     }


   rotAngle[0]=-1.*gammaList[0];
   rotAngle[1]=-1.*betaList[0];
   rotAngle[2]=-1.*alphaList[0];

   if( error = rotateDomain(rotAngle,pdbInfo,pdbInfo->minResID,pdbInfo->maxResID) )
     {
	FPR( stderr, "matrixt: error in rotateDomain." );
	goto shutdown;
     }


   if( *pdbRotName )
     {
        (void) writePDB( pdbRotName, pdbInfo, PDB_SELECT_NULL, remark );
     }


   if( error = restorePDBCoords( pdbInfo ) )
     {
         FPR( stderr, "matrix: error in storePDBCoords.\n" );
         goto shutdown;
     }



/***/
/* Exit point:
/***/

shutdown:

    return (error); 
}
/* ========================================================================= */


/* ========================================================================= */
int selectSurfDom( struct PDBInfo *pdbInfo, int resN, int resC,
	 int *indICount, int *indI, int *indJCount, int *indJ )
{
   int   i,error,selectCount;
   int   dom1Count,dom2Count;


/***/
/* Initialize:
/***/

   i           = 0;
   error       = 0;
   selectCount = 0;
   dom1Count   = 0;
   dom2Count   = 0;


/***/
/* Select surface accessible of rotated domain:
/***/

   selectCount = pdbSetSelection( pdbInfo, 0 );
   selectCount = pdbSelect1( pdbInfo, resN, resC, "C", 1, 3 );
   selectCount = pdbSelect1( pdbInfo, resN, resC, "CA", 1, 3 );
   selectCount = pdbSelect1( pdbInfo, resN, resC, "N", 1, 3 );
   if( selectCount==0 )
     {
	FPR( stderr, "Error in chkDomClash: zero atoms selected." );
	error = 1;
	goto shutdown;
     }
   if ( error = callNSC( pdbInfo ) )
      {
        goto shutdown;
      }


/***/
/* Store selected coordinates of domain 1 in coordinate structure:
/***/

  dom1Count = 0; 
  for( i = 0; i < pdbInfo->atomCount; i++ )
       {
       if( pdbInfo->atomList[i].selected )
          {
           indI[dom1Count] = i; 
           dom1Count++;
          }
       }
  *indICount = dom1Count;



/***/
/* Select surface accessible atoms of other part of molecule:
/***/

   selectCount = pdbSetSelection( pdbInfo, 0 );
   selectCount = pdbSelect1( pdbInfo, pdbInfo->minResID, pdbInfo->maxResID, "C", 1, 3 );
   selectCount = pdbSelect1( pdbInfo, pdbInfo->minResID, pdbInfo->maxResID, "N", 1, 3 );
   selectCount = pdbSelect1( pdbInfo, pdbInfo->minResID, pdbInfo->maxResID, "CA", 1, 3 );
   selectCount = pdbSelect1( pdbInfo, resN, resC, "*", 0, 3 );
   if( selectCount==0 )
     {
	FPR( stderr, "Error in chkDomClash: zero atoms selected." );
	error = 1;
	goto shutdown;
     }

   if ( error = callNSC( pdbInfo ) )
      {
        goto shutdown;
      }


/***/
/* Store selected coordinates of domain 1 in coordinate structure:
/***/

  dom2Count = 0; 
  for( i = 0; i < pdbInfo->atomCount; i++ )
       {
       if( pdbInfo->atomList[i].selected )
          {
           indJ[dom2Count] = i; 
           dom2Count++;
          }
       }
  *indJCount = dom2Count;



/***/
/* Exit point:
/***/

shutdown:

    return (error); 
}
/* ========================================================================= */


/* ============================================================ */
float calcNorm( float x,  float y,  float z )
{
    float norm=0;

    norm = sqrt( x*x + y*y + z*z );

    return ( norm );
}
/* ============================================================ */


/* ============================================================ */
int ssiaFree( struct PalesInfo *palesInfo, char *pdbFreeName,
                char *outAName, struct DTable *dTab )
{
   char    result[6000];
   int     i, error, density;
   int     atomCount,atomSelCount,dotCount;
   float   center_of_mass[3];
   float   maxRadius, minRadius, distance;
   double  *xI,*yI,*zI;
   double  *dots,*ddots;
   double  saupe[5];
   double  radius,gridSpacing,wvAmount,molarRatio;



/***/
/* Initialize:
/***/

    error        = 0;
    i            = 0;
    density      = 0;
    dotCount     = 0;
    distance     = 0;
    maxRadius    = 0.0;
    minRadius    = 100000.0;
    atomCount    = 0;
    atomSelCount = 0;


    for( i=0; i<3; i++ )     center_of_mass[i]=0;


    *result   = '\0';
   
    xI           = (double *)NULL;
    yI           = (double *)NULL;
    zI           = (double *)NULL;

    dots         = (double *)NULL;
    ddots        = (double *)NULL;


/***/
/* Allocate and extract coord lists, etc:
/***/

    xI    = dFltAlloc( "ssiaFree", MAX_ATOM_NO+1 );
    yI    = dFltAlloc( "ssiaFree", MAX_ATOM_NO+1 );
    zI    = dFltAlloc( "ssiaFree", MAX_ATOM_NO+1 );
    ddots = dFltAlloc( "ssiaFree", (MAX_ENSEMBLE+1)*3 );

    if (!xI || !yI || !zI || !ddots)
       {
        FPR( stderr, "Error allocating coord buffers.\n" );
        error = 1;
        goto shutdown;
       }



   refine_bicelle( palesInfo );


/***/
/* Read PDB file in free format:
/***/

   if( error = readPDBfree( pdbFreeName,&atomCount,xI,yI,zI ) )
     {
        FPR( stderr, "Matrix: Error reading PDB file in free format.\n" );
        goto shutdown;
     }    

   FPR( stderr,
        "%s, %d Atoms, Free format.\n", pdbFreeName, atomCount );



/***/
/* Setup grid of orientations:
/***/

   density = palesInfo->densityNSC;
   if( UNSP(density,&dots,&dotCount) )
     {
        FPR( stderr, "Matrix: unsp error occured.\n" );
	error = 1;
        goto shutdown;
     }
   for( i=0; i<dotCount; i++ )
      {
        ddots[i*3]   = dots[i*3];
        ddots[i*3+1] = dots[i*3+1];
        ddots[i*3+2] = dots[i*3+2];
      }


/***/
/* Select surface accessible atoms:
/***/

   if( error = surfSel(atomCount,xI,yI,zI,&atomSelCount) )
     {
        FPR( stderr, "Matrix: Error reading PDB file in free format.\n" );
        goto shutdown;
     }    

   if( atomSelCount==0 )
     {
        FPR( stderr, "PALES: error - zero atoms remain after surface selection!\n" );
        error = 1;
        goto shutdown;
     }    

   FPR( stderr, "REMARK %d atoms selected for simulation.\n", atomSelCount);


/***/
/* Shift protein to geometric center
/* Determine maximum and minimum radius of protein
/***/

   /*** calculate geometric center of positioning helix ***/
   for( i=0; i<3; i++ )     center_of_mass[i]=0;
   for( i=1; i<=atomSelCount; i++ )
      {
                center_of_mass[0] += xI[i];
                center_of_mass[1] += yI[i];
                center_of_mass[2] += zI[i];
      }

   center_of_mass[0] /= atomSelCount;
   center_of_mass[1] /= atomSelCount;
   center_of_mass[2] /= atomSelCount;


  /*** shift selected part of protein to geometric center ***/
  for( i=1; i<=atomSelCount; i++ )
     {
         xI[i]-= center_of_mass[0];
         yI[i]-= center_of_mass[1];
         zI[i]-= center_of_mass[2];
     }


  /*** determine extrema of length of protein ***/
  minRadius = 100000.0;
  maxRadius = 0.0;
  for( i=1; i<=atomSelCount; i++ )
     {
        distance = xI[i]*xI[i]+yI[i]*yI[i]+zI[i]*zI[i];
        if( distance > maxRadius ) maxRadius = distance;
        if( distance < minRadius ) minRadius = distance;
     }
  minRadius = sqrt( minRadius );
  maxRadius = sqrt( maxRadius );


/***/
/* Simulate sterically induced alignment:
/***/

   radius      = palesInfo->radiusAtom;
   gridSpacing = palesInfo->gridSpacing;
   wvAmount    = palesInfo->concentration;
   molarRatio  = palesInfo->molarRatio;

   if( error = stericAlignment( palesInfo, palesInfo->orderDiscRim,palesInfo->lcType,palesInfo->rSingleCell,palesInfo->radiusModel,
              radius,gridSpacing,wvAmount,molarRatio,
	      atomSelCount,xI,yI,zI,maxRadius+palesInfo->radiusAtom,minRadius+palesInfo->radiusAtom,
              palesInfo->digPsiNSC,dotCount,ddots,saupe) )
     {
	FPR( stderr, "simulateAlignment: error in stericAlignment." );
	goto shutdown;
     } 

/*
  ssia( palesInfo,atomSelCount,xI,yI,zI,
        maxRadius+palesInfo->radiusAtom,minRadius+palesInfo->radiusAtom,
        dotCount,ddots,saupe );
*/

  /* scale Saupe matrix according to bicelle order */
  for( i=0; i<5; i++ )    saupe[i] *= palesInfo->lcOrder; 


/***/
/* Selection of sign according to Frank's convention:
/***/

   for( i = 0; i < 5; i++ )
      {
         saupe[i] = (-1.0)*saupe[i];
      }



/***/
/* Write results:
/***/

   dTab->dTensor.saupe[0] = saupe[0]; 
   dTab->dTensor.saupe[1] = saupe[1]; 
   dTab->dTensor.saupe[2] = saupe[2]; 
   dTab->dTensor.saupe[3] = saupe[3]; 
   dTab->dTensor.saupe[4] = saupe[4]; 

/*
   dTab->dTensor.eulerCount = rot2euler( eigVec, 
                                     dTab->dTensor.psiList, 
                                     dTab->dTensor.thetaList, 
                                     dTab->dTensor.phiList );
*/
 

   (void) writeDCmod( outAName, "w", dTab, palesInfo, 
               result );



/***/
/* Exit point:
/***/

shutdown:
 
    (void) deAlloc( "ssiaFree", xI,    sizeof(double)*(MAX_ATOM_NO+1) );
    (void) deAlloc( "ssiaFree", yI,    sizeof(double)*(MAX_ATOM_NO+1) );
    (void) deAlloc( "ssiaFree", zI,    sizeof(double)*(MAX_ATOM_NO+1) );

    (void) deAlloc( "ssiaFree", ddots, sizeof(double)*(MAX_ENSEMBLE+1)*3 );

    return( error );

}
/* ============================================================ */


/* ============================================================ */
int shiftMol2Domain( struct PDBInfo *pdbInfo, int resN, int resC )
{
  int     i,atomCount,selectCount,error;
  float   *xI, *yI, *zI;
  float   center_of_mass[3];


/***/
/* Initialize:
/***/

    error       = 0;
    i           = 0;
    selectCount = 0;

    xI        = (float *)NULL;
    yI        = (float *)NULL;
    zI        = (float *)NULL;


    atomCount = pdbInfo->atomCount;


/***/
/* Allocate memory:
/***/

    xI = fltAlloc( "shiftMol2Domain", atomCount );
    yI = fltAlloc( "shiftMol2Domain", atomCount );
    zI = fltAlloc( "shiftMol2Domain", atomCount );
    if (!xI || !yI || !zI)
       {
        FPR( stderr, "Error in shiftMol2Domain: allocating coord buffers.\n" );
        error = 1;
        goto shutdown;
       }



/***/
/* Calculate geometric center of selected part of molecule:
/***/

 for( i = 0; i < pdbInfo->atomCount; i++ )
     {
        if( pdbInfo->atomList[i].resID >= resN )
          {
           if( pdbInfo->atomList[i].resID <= resC ) 
             {
               xI[selectCount]    = pdbInfo->atomList[i].x;
               yI[selectCount]    = pdbInfo->atomList[i].y;
               zI[selectCount]    = pdbInfo->atomList[i].z;
               selectCount++;
             }
           else
                break;
           }
     }
  if( selectCount==0 )
    {
	FPR( stderr, "Error in positionDomain: zero atoms selected." );
	error = 1;
	goto shutdown;
    }


   /*** calculate geometric center ***/
   for( i=0; i<3; i++ )     center_of_mass[i] = 0;
   for( i=0; i<selectCount; i++ )
      {
                center_of_mass[0] += xI[i];
                center_of_mass[1] += yI[i];
                center_of_mass[2] += zI[i];
      }

   center_of_mass[0] /= selectCount;
   center_of_mass[1] /= selectCount;
   center_of_mass[2] /= selectCount;



/***/
/* Shift whole molecule to geometric center of selected part:
/***/

  for( i = 0; i < pdbInfo->atomCount; i++ )
     {
        pdbInfo->atomList[i].x -= center_of_mass[0];
        pdbInfo->atomList[i].y -= center_of_mass[1];
        pdbInfo->atomList[i].z -= center_of_mass[2];
     }



/***/
/* Exit point:
/***/

shutdown:
 
    (void) deAlloc( "shiftMol2Domain", xI,      sizeof(double)*atomCount );
    (void) deAlloc( "shiftMol2Domain", yI,      sizeof(double)*atomCount );
    (void) deAlloc( "shiftMol2Domain", zI,      sizeof(double)*atomCount );

    return( error );
}
/* ========================================================================= */


/* ========================================================================= */
int  mapSimDC( struct PDBInfo *pdbInfo, struct DTable *dTab, 
           struct PalesInfo *palesInfo, 
           char *outDName, char *pdbRotName, int verbose,
           char *outMapName, int mapCount, char *remark)
{
   char    ctemp[6000], result[6000];
   int     i,j,n,error,fi,cycle1,cycle3;
   float   rSum, pSum;
   float   *dcSum , *dcSqrSum, *dcAve, *dcVar, *dcStdDev;
   float   *dVal, *dSigma;
   int     *countFail;
   double  szzSum, szzSqrSum, szzAve, szzVar, stdDev;
   double  fac, kappa;
   double  szz, sxx_yy, sxy, sxz, syz, sxx, syy;
   double  gMagTemp, gMagSum, gMagSqrSum, gMagAve, gMagVar, gMagStdDev;
   double  ran, lambda, phi, mapx, mapy;
   double  rotAngle[3];
   double  eigVal[3], eigVec[3][3];
   double  gMag, da, dr, saupe[5], irrSaupe[5];
   double  saupeSum[5], saupeSqrSum[5], saupeAve[5], saupeVar[5], saupeStdDev[5];
   static  double alphaList[NO_LIST], betaList[NO_LIST], gammaList[NO_LIST];
   FILE    *outUnit;



/***/
/* Initialize:
/***/

   i         = 0;
   j         = 0;
   n         = 0;
   error     = 0;
   outUnit   = stdout;

   cycle1    = 0;
   cycle3    = 0;

   rSum      = 0;
   pSum      = 0;
   gMag      = 0;
   da        = 0;
   dr        = 0;

   szzSum    = 0;
   szzSqrSum = 0;   
   szzAve    = 0;
   szzVar    = 0;
   stdDev    = 0;


   kappa  = 21523;
   fac    = (2.0/3.0) * (4.0/5.0*PI) * kappa*kappa;

   szz    = 0;
   szz    = 0;
   sxx_yy = 0;
   sxy    = 0;
   sxz    = 0;
   syz    = 0;
   sxx    = 0;
   syy    = 0;

   gMagTemp   = 0;
   gMagSum    = 0;
   gMagSqrSum = 0;
   gMagAve    = 0;
   gMagVar    = 0; 
   gMagStdDev = 0;


   ran       = 0;
   lambda    = 0;
   phi       = 0;
   mapx      = 0;
   mapy      = 0;

   fi        = 1;


   *ctemp    = '\0';
   *result   = '\0';


   dcSum     = (float *)NULL;
   dcSqrSum  = (float *)NULL;
   dcAve     = (float *)NULL;
   dcVar     = (float *)NULL;
   dcStdDev  = (float *)NULL;

   dVal      = (float *)NULL;
   dSigma    = (float *)NULL;
   countFail = (int *)NULL;
 



   for( i=0; i<5; i++ )
     {
          saupe[i]       = 0.0; 
          irrSaupe[i]    = 0.0; 

          saupeSum[i]    = 0.0; 
          saupeSqrSum[i] = 0.0; 
          saupeAve[i]    = 0.0; 
          saupeVar[i]    = 0.0; 
          saupeStdDev[i] = 0.0; 
     }


   if (!dTab)            return( 0 );
   if (!dTab->dCount)    return( 0 );

   if (!dTab->dObsList)  return( 0 );
   if (!dTab->dCalcList) return( 0 );




/***/
/* Allocate array for average couplings:
/***/

  dcSum     = fltAlloc( "matrix", dTab->dCount );
  dcSqrSum  = fltAlloc( "matrix", dTab->dCount );
  dcAve     = fltAlloc( "matrix", dTab->dCount );
  dcVar     = fltAlloc( "matrix", dTab->dCount );
  dcStdDev  = fltAlloc( "matrix", dTab->dCount );

  dVal      = fltAlloc( "matrix", dTab->dCount );
  dSigma    = fltAlloc( "matrix", dTab->dCount );

  countFail = intAlloc( "matrix", dTab->dCount );



  if (!dcAve || !dcSum || !dcSqrSum || !dcVar || !dcStdDev)
     {
        FPR( stderr, "Error allocating coord buffers.\n" );
        error = 1;
        goto shutdown;
     }


  if (!dVal || !dSigma || !countFail)
     {
        FPR( stderr, "Error allocating coord buffers.\n" );
        error = 1;
        goto shutdown;
     }


  for( i = 0; i < dTab->dCount; i++ )
     {
        dcSum[i]     = 0;
        dcSqrSum[i]  = 0;
        dcAve[i]     = 0;
        dcVar[i]     = 0;
        dcStdDev[i]  = 0;

        dVal[i]      = 0;
        dSigma[i]    = 0;

        countFail[i] = 0;
     }



/***/
/* Open the output files:
/***/

   outUnit = stdout;

   if ( *outMapName )
      {
        if (!(outUnit = fopen( outMapName, "w" )))
           {
            FPR( stderr, "PALES Error opening log file for error space mapping %s\n", outMapName );
            outUnit = stdout;
            error = 1;
            goto shutdown;
           }
      }


 
/***/
/* Store observed DCs and their errors:
/***/

   for( i = 0; i < dTab->dCount; i++ )
      {
         if (!dTab->dCalcList[i].valid) continue;
         if (!dTab->dCalcList[i].atomI[0]) continue;
         if (!dTab->dCalcList[i].atomJ[0]) continue;

         dVal[i]   = dTab->dObsList[i].d;
         dSigma[i] = dTab->dObsList[i].sigma;
      }       


/***/
/* Map error space:
/***/

   for( i=0; i<mapCount; i++ )
      {
       	 cycle1 = cycle1 + 1;
         fi     =  1;


/***/
/* genarating random numbers with a gaussian distribution, with a mean of 0 and a variance of 1  */
/* this way, the variance of our distribution is error/2 , the std is the error                 */
/* See "Numerical Recepies in C", Ch 7                                                          */
/***/

         for( j = 0; j < dTab->dCount; j++ )
            {
               if (!dTab->dCalcList[j].valid) continue;
               if (!dTab->dCalcList[j].atomI[0]) continue;
               if (!dTab->dCalcList[j].atomJ[0]) continue;

               ran = gauss();

               dTab->dObsList[j].d     = dVal[j] + (ran * dSigma[j]/2.0);

            }       

/***/
/* Best-fit dipolar couplings:
/***/

         (void) fitD( pdbInfo, dTab, verbose );


/***/
/*  Checking if values are within the experimental range:
/***/


         for( j = 0; j < dTab->dCount; j++ )
           {
              if (!dTab->dCalcList[j].valid) continue;
              if (!dTab->dCalcList[j].atomI[0]) continue;
              if (!dTab->dCalcList[j].atomJ[0]) continue;

/*
	      if( (dTab->dCalcList[j].d > dVal[j] + dSigma[j]) || (dTab->dCalcList[j].d < dVal[j] - dSigma[j]) )
                {
	           fi = 0;
	           countFail[j]++;
                }
*/

	      fi = 0;
	      if( fabs(dTab->dCalcList[j].d - dVal[j]) < dSigma[j] && fabs(dTab->dCalcList[j].d - dVal[j]) >= 0 ) 
                {
	           fi = 1;
                }
              if( fi == 0 )
	           countFail[j]++;
           }


/***/
/*  fi == 1 only if the solutions are within experimental error  */
/***/


         if( fi == 1 )
           { 
  	      cycle3++;

              for( j = 0; j < 5; j++ )
                 {
                   saupe[j] = dTab->dTensor.saupe[j];
                 }


              if( analyzeSaupe(saupe, irrSaupe, &gMag, eigVec, eigVal, &da, &dr,
                 alphaList, betaList, gammaList, palesInfo->eulerFlag) )
                {
                   FPR( stderr, "Matrix: error in analyzeSaupe.\n");
                   error = 1;
                   goto shutdown;
                }


              for( j=0; j<5; j++ )
                {
                    saupeSum[j]    += saupe[j]; 
                    saupeSqrSum[j] += saupe[j]*saupe[j]; 
                }

	      szzSum     += fabs(eigVal[0]);
	      szzSqrSum  += eigVal[0]*eigVal[0];


              szz    = saupe[0];
              sxx_yy = saupe[1];
              sxy    = saupe[2];
              sxz    = saupe[3];
              syz    = saupe[4];
              sxx    = -0.5*(szz - sxx_yy);
              syy    = -0.5*(szz + sxx_yy);

              gMagTemp    =  sqrt( fac *(sxx*sxx + syy*syy + szz*szz + 2*sxy*sxy + 2*sxz*sxz + 2*syz*syz ));
              gMagSum    += gMagTemp;
              gMagSqrSum += gMagTemp*gMagTemp;


	      for( j = 0; j < dTab->dCount; j++ )
      		{
                   if (!dTab->dCalcList[j].valid) continue;
                   if (!dTab->dCalcList[j].atomI[0]) continue;
                   if (!dTab->dCalcList[j].atomJ[0]) continue;

                   dcSum[j]    += dTab->dCalcList[j].d;
                   dcSqrSum[j] += dTab->dCalcList[j].d*dTab->dCalcList[j].d;
                }


/***/
/*  Output of mapping:
/***/

	      phi    = PI/2-acos(eigVec[2][0]);
	      lambda = atan2(eigVec[1][0], eigVec[0][0]);
	      mapx   = lambda*cos(phi);
	      mapy   = phi;
	      FPR( outUnit, "  %f       %f  ", mapx, mapy); 
	       		
	      phi    = PI/2-acos(eigVec[2][1]);
	      lambda = atan2(eigVec[1][1], eigVec[0][1]);
	      mapx   = lambda*cos(phi);
	      mapy   = phi;
	      FPR( outUnit, "  %f       %f  ", mapx, mapy); 
				 
	      phi    = PI/2-acos(eigVec[2][2]);
	      lambda = atan2(eigVec[1][2], eigVec[0][2]);
	      mapx   = lambda*cos(phi);
	      mapy   = phi;
	      FPR( outUnit, "  %f       %f\n", mapx, mapy); 
	      
	      phi    = PI/2-acos(eigVec[2][0]);
	      lambda = atan2(eigVec[1][0], eigVec[0][0]);
	      mapx   = lambda*cos(phi);
	      mapy   = phi;
	      if (lambda>0) mapx=(lambda-PI)*cos(-phi);
	      if (lambda<0) mapx=(lambda+PI)*cos(-phi);
	      mapy   = -phi;
	      FPR( outUnit, " %f	%f  ", mapx, mapy);  
	       		
	      phi    = PI/2-acos(eigVec[2][1]);
	      lambda = atan2(eigVec[1][1], eigVec[0][1]);
	      mapx   = lambda*cos(phi);
	      mapy   = phi;
	      if (lambda>0) mapx=(lambda-PI)*cos(-phi);
	      if (lambda<0) mapx=(lambda+PI)*cos(-phi);
	      mapy  = -phi;
	      FPR( outUnit, " %f	%f  ", mapx, mapy);  
				 
	      phi    = PI/2-acos(eigVec[2][2]);
	      lambda = atan2(eigVec[1][2], eigVec[0][2]);
	      mapx   = lambda*cos(phi);
	      mapy   = phi;
	      if (lambda>0) mapx=(lambda-PI)*cos(-phi);
	      if (lambda<0) mapx=(lambda+PI)*cos(-phi);
	      mapy = -phi;
	      FPR( outUnit, " %f	%f\n", mapx, mapy); 


	   } /*  end of error space search  */

      } /* end of mapCount loop */
 

   FPR( stdout, "\nREMARK Number of steps for error statistic %d.\n", cycle1 );



/***/	
/* Calculate & store average values:
/***/	

   if( cycle3 != 0 )
     {
        szzAve  = szzSum    / cycle3;
        szzVar  = szzSqrSum / cycle3 - szzAve*szzAve;
        stdDev  = sqrt( szzVar );

        gMagAve     = gMagSum    / cycle3;
        gMagVar     = gMagSqrSum / cycle3 - gMagAve*gMagAve;
        gMagStdDev  = sqrt( gMagVar );
     }



/***/	
/* Calculate & store average values:
/***/	

   if( cycle3 != 0 )
     {

       for( i = 0; i < 5; i++ )
         {
            saupeAve[i]    = saupeSum[i]    / cycle3;
            saupeVar[i]    = saupeSqrSum[i] / cycle3 - saupeAve[i]*saupeAve[i];
            saupeStdDev[i] = sqrt( saupeVar[i] );

            dTab->dTensor.saupe[i] = saupeAve[i];
         }


       for( i = 0; i < dTab->dCount; i++ )
         {
            if (!dTab->dCalcList[i].valid) continue;
            if (!dTab->dCalcList[i].atomI[0]) continue;
            if (!dTab->dCalcList[i].atomJ[0]) continue;

            dcAve[i]    = dcSum[i]    / cycle3;
            dcVar[i]    = dcSqrSum[i] / cycle3 - dcAve[i]*dcAve[i];
            dcStdDev[i] = sqrt( dcVar[i] );

            dTab->dObsList[i].d     = dVal[i];
            dTab->dCalcList[i].d    = dcAve[i];
            dTab->dObsList[i].sigma = dcStdDev[i];
            dTab->dObsList[i].w     = countFail[i];
         }





/***/
/* Analyze average back-calculated DCs:
/***/

        SPR( result, "\nREMARK Mapping Statistics.\n\n" );

        *ctemp    = '\0';
        SPR( ctemp, "DATA STATISTICS MAPPING NUMBER OF SUCCESFUL STEPS   %d\n", cycle3 );
        (void) strcat( result, ctemp );

        *ctemp    = '\0';
        SPR( ctemp, "DATA STATISTICS MAPPING G_MAG (AVERAGE, STD_DEV)   % .4e % .4e\n", gMagAve, gMagStdDev );   
        (void) strcat( result, ctemp );

        *ctemp    = '\0';
        SPR( ctemp, "DATA STATISTICS MAPPING SZZ (AVERAGE, STD_DEV)     % .4e % .4e\n", szzAve, stdDev );   
        (void) strcat( result, ctemp );

        *ctemp    = '\0';
        SPR( ctemp, "DATA STATISTICS MAPPING SAUPE AVERAGE    % .4e  % .4e  % .4e  % .4e % .4e\n", 
        saupeAve[0], saupeAve[1], saupeAve[2], saupeAve[3], saupeAve[4]);   
        (void) strcat( result, ctemp );

        *ctemp    = '\0';
        SPR( ctemp, "DATA STATISTICS MAPPING SAUPE STD_DEV    % .4e  % .4e  % .4e  % .4e % .4e\n", 
        saupeStdDev[0], saupeStdDev[1], saupeStdDev[2], saupeStdDev[3], saupeStdDev[4]);   
        (void) strcat( result, ctemp );

        *ctemp    = '\0';
        SPR( ctemp, "\n\n" );
        (void) strcat( result, ctemp );


/***/
/* Write results:
/***/

        (void) writeDCmod( outDName, "w", dTab, palesInfo, 
                    result );



/***/
/* Write molecule rotated such that alignment tensor coincide with laboratory frame:
/***/

       if( error = storePDBCoords( pdbInfo ) )
         {
            FPR( stderr, "matrix: error in storePDBCoords.\n" );
            goto shutdown;
         }


       rotAngle[0] = -1.*gammaList[0];
       rotAngle[1] = -1.*betaList[0];
       rotAngle[2] = -1.*alphaList[0];

       if( error = rotateDomain(rotAngle,pdbInfo,pdbInfo->minResID,pdbInfo->maxResID) )
         {
	    FPR( stderr, "matrixt: error in rotateDomain." );
	    goto shutdown;
         }


       if( *pdbRotName )
         {
            (void) writePDB( pdbRotName, pdbInfo, PDB_SELECT_NULL, remark );
         }


       if( error = restorePDBCoords( pdbInfo ) )
         {
         FPR( stderr, "matrix: error in storePDBCoords.\n" );
         goto shutdown;
         }

     } /* end of if cycle3!=0 loop */


/***/
/* Exit point:
/***/

shutdown:

    if (outUnit != stdout) (void) fclose( outUnit );

    (void) deAlloc( "matrix", dcSum,     sizeof(float)*dTab->dCount );
    (void) deAlloc( "matrix", dcSqrSum,  sizeof(float)*dTab->dCount );
    (void) deAlloc( "matrix", dcAve,     sizeof(float)*dTab->dCount );
    (void) deAlloc( "matrix", dcVar,     sizeof(float)*dTab->dCount );
    (void) deAlloc( "matrix", dcStdDev,  sizeof(float)*dTab->dCount );

    (void) deAlloc( "matrix", dVal,      sizeof(float)*dTab->dCount );
    (void) deAlloc( "matrix", dSigma,    sizeof(float)*dTab->dCount );

    (void) deAlloc( "matrix", countFail, sizeof(int)*dTab->dCount );


    return (error); 
}
/* ========================================================================= */


/* ========================================================================= */
int compareSaupe( struct PalesInfo *palesInfo,
                double userSaupe1[5], double userSaupe2[5],
		char *outName, struct DTable *dTab, char *remark)
{
   char    ctemp[6000], result[6000];
   int     i,error;
   float   dcos[3], dcos5d, angle;
   double  gMag, gMag1, gMag2, R;
   double  saupe[5], irrSaupe[5], irrSaupe1[5], irrSaupe2[5];
   double  da, dr;
   double  da1, dr1, da2, dr2;
   double  eigVal[3], eigVec[3][3];
   double  eigVal1[3], eigVec1[3][3];
   double  eigVal2[3], eigVec2[3][3];
   static  double alphaList[NO_LIST], betaList[NO_LIST], gammaList[NO_LIST];
   static  double alphaList1[NO_LIST], betaList1[NO_LIST], gammaList1[NO_LIST];
   static  double alphaList2[NO_LIST], betaList2[NO_LIST], gammaList2[NO_LIST];



/***/
/* Initialize:
/***/

   i         = 0;
   error     = 0;

   gMag      = 0;
   da        = 0;
   dr        = 0;
   gMag1     = 0;
   gMag2     = 0;
   da1       = 0;
   dr1       = 0;
   da2       = 0;
   dr2       = 0;
   R         = 0;
   dcos5d    = 0.0;
   angle     = 0.0;

   *ctemp    = '\0';
   *result   = '\0';



/***/
/* Apply operation to two order matrices:
/***/

    if( addSaupe(palesInfo->opCode, palesInfo->c1, palesInfo->c2,
                  userSaupe1, userSaupe2, saupe) )
      {
        FPR( stderr, "pales: error in addSaupe.\n");
        error = 1;
        goto shutdown;
      }



/***/
/* Analyze individual order matrices:
/***/

    if( analyzeSaupe(userSaupe1, irrSaupe1, &gMag1, eigVec1, eigVal1, &da1, &dr1,
        alphaList1, betaList1, gammaList1, palesInfo->eulerFlag) )
      {
        FPR( stderr, "Matrix: error in analyzeSaupe.\n");
        error = 1;
        goto shutdown;
       }

    if( analyzeSaupe(userSaupe2, irrSaupe2, &gMag2, eigVec2, eigVal2, &da2, &dr2,
        alphaList2, betaList2, gammaList2, palesInfo->eulerFlag) )
      {
        FPR( stderr, "Matrix: error in analyzeSaupe.\n");
        error = 1;
        goto shutdown;
       }

    if( analyzeSaupe(saupe, irrSaupe, &gMag, eigVec, eigVal, &da, &dr,
        alphaList, betaList, gammaList, palesInfo->eulerFlag) )
      {
        FPR( stderr, "Matrix: error in analyzeSaupe.\n");
        error = 1;
        goto shutdown;
      }



/***/
/* Compare order matrices:
/***/

   (void) calcDcos( eigVec1, eigVec2, dcos );


   if(  irrCollinarity( userSaupe1, userSaupe2, &dcos5d, &angle ) )
      {
        FPR( stderr, "pales: error in compareSaupe.\n");
        error = 1;
        goto shutdown;
      }



/***/
/* Output:
/***/

   SPR( result,
         "\nREMARK Input order matrices.\n\n" );

   *ctemp    = '\0';
   SPR( ctemp,
            "DATA SAUPE 1      % .4e % .4e % .4e % .4e % .4e\n",
            userSaupe1[0], 
            userSaupe1[1], 
            userSaupe1[2], 
            userSaupe1[3], 
            userSaupe1[4] );
   (void) strcat( result, ctemp );

   *ctemp    = '\0';
   SPR( ctemp,
            "DATA SAUPE 2      % .4e % .4e % .4e % .4e % .4e\n",
            userSaupe2[0], 
            userSaupe2[1], 
            userSaupe2[2], 
            userSaupe2[3], 
            userSaupe2[4] );
   (void) strcat( result, ctemp );

   *ctemp    = '\0';
   SPR( ctemp,
            "DATA SAUPE RESULT % .4e % .4e % .4e % .4e % .4e\n\n",
            saupe[0], 
            saupe[1], 
            saupe[2], 
            saupe[3], 
            saupe[4] );
   (void) strcat( result, ctemp );


   *ctemp    = '\0';
   SPR( ctemp,
         "\nREMARK Relative Orientation of alignment tensors.\n\n" );
   (void) strcat( result, ctemp );


   *ctemp    = '\0';
   SPR( ctemp,
         "DATA ANGLE_3D_AXES (X/Y/Z)  %.2f %.2f %.2f\n",
         DEG(dcos[0]), DEG(dcos[1]), DEG(dcos[2]) );
   (void) strcat( result, ctemp );

   *ctemp    = '\0';
   SPR( ctemp,
         "\nDATA ANGLE_5D_SPACE  %.2f\n", angle );
   (void) strcat( result, ctemp );

   *ctemp    = '\0';
   SPR( ctemp,
         "DATA COLL_5D         %.5f\n\n\n", dcos5d );
   (void) strcat( result, ctemp );



   dTab->dTensor.saupe[0] = saupe[0]; 
   dTab->dTensor.saupe[1] = saupe[1]; 
   dTab->dTensor.saupe[2] = saupe[2]; 
   dTab->dTensor.saupe[3] = saupe[3]; 
   dTab->dTensor.saupe[4] = saupe[4]; 


   (void) writeDCmod( outName, "w", dTab, palesInfo, 
                 result);



/***/
/* Exit point:
/***/

shutdown:

    return (error); 
}
/* ========================================================================= */

/* ========================================================================= */
int addSaupe( int opCode, float c1, float c2,
                double userSaupe1[5], double userSaupe2[5],
                double saupe[5] )
{
   int     i;
   double  gMag1, gMag2;
   double  irrSaupe1[5], irrSaupe2[5];


/***/
/* Initialize:
/***/

   i         = 0;
   gMag1     = 0.0;
   gMag2     = 0.0;

   for( i = 0; i < 5; i++ )
     {
        saupe[i] = 0;
        irrSaupe1[i] = 0;
        irrSaupe2[i] = 0;
     }



/***/
/* Apply operation to two order matrices:
/***/


   switch( opCode )
      {
       case OP_ADD:
          for( i = 0; i < 5; i++ )
            {
               saupe[i] = c1 * userSaupe1[i] + c2 * userSaupe2[i];
            }
          break;

       case OP_SUB:
          for( i = 0; i < 5; i++ )
            {
               saupe[i] = c1 * userSaupe1[i] - c2 * userSaupe2[i];
            }
          break;

       case OP_MUL:
          for( i = 0; i < 5; i++ )
            {
               saupe[i] = c1 * userSaupe1[i] * c2 * userSaupe2[i];
            }
          break;

       case OP_DIV:
          for( i = 0; i < 5; i++ )
            {
               saupe[i] = c1 * userSaupe1[i] / (c2 * userSaupe2[i]);
            }
          break;

       case OP_MAX:
          irreducRepres( userSaupe1, &gMag1, irrSaupe1 );
          irreducRepres( userSaupe2, &gMag2, irrSaupe2 );

          if( gMag1 >= gMag2 )
            {
               for( i = 0; i < 5; i++ )
                 {
                    saupe[i] = c1 * userSaupe1[i];
                 }
            }
          else
            {
               for( i = 0; i < 5; i++ )
                 {
                    saupe[i] = c2 * userSaupe2[i];
                 }
            }
          break;
 
       default:
          break;
      }


/***/
/* Exit point:
/***/

    return( 0 );
}
/* ========================================================================= */


/* ========================================================================= */
int analyzeDC( char *inD1Name, char *inD2Name,
               struct PalesInfo *palesInfo, struct DcInfo *dcInfo,
	       char *atomNameS1, char *atomNameS2, char *outName)
{
   char            result[6000], ctemp[6000];
   int             i,error,validCountD;
   int             itemp;
   struct DTable   *dTabAn1, *dTabAn2;


/***/
/* Initialize:
/***/

   i           = 0;
   error       = 0;
   validCountD = 0;

   *ctemp    = '\0';
   *result   = '\0';

   dTabAn1     = (struct DTable *)NULL;
   dTabAn2     = (struct DTable *)NULL;



/***/
/* Allocate PDB and DC structures.
/* Null-initialize them.
/***/

    dTabAn1    = (struct DTable *)voidAlloc( "pales",  sizeof(struct DTable) );
    dTabAn2    = (struct DTable *)voidAlloc( "pales",  sizeof(struct DTable) );

    if (!dTabAn1 || !dTabAn1)
       {
        FPR( stderr, "Error allocating DC structures.\n" );
        error = 1;
        goto shutdown;
       }

    (void) dTabNull( dTabAn1 );
    (void) dTabNull( dTabAn2 );



/***/
/* Read the input files:
/***/

   if (error = readDC( inD1Name, dTabAn1 ))
      {
        FPR( stderr, "Error reading DC input %s\n", inD1Name );
        goto shutdown;
      }

   if (error = readDC( inD2Name, dTabAn2 ))
      {
        FPR( stderr, "Error reading DC input %s\n", inD2Name );
        goto shutdown;
      }

   FPR( stderr,
          "%s, %d Couplings, Residues %d to %d.\n",
           inD1Name, dTabAn1->dCount, dTabAn1->dMinRes, dTabAn1->dMaxRes );

   FPR( stderr,
          "%s, %d Couplings, Residues %d to %d.\n",
          inD2Name, dTabAn2->dCount, dTabAn2->dMinRes, dTabAn2->dMaxRes );


   if (dcInfo->segFirst < -99999)   dcInfo->segFirst = dTabAn1->dMinRes;
   if (dcInfo->segLast  < -99999)   dcInfo->segLast  = dTabAn1->dMaxRes;

   if (dcInfo->segFirst > dcInfo->segLast) {ISWAP( dcInfo->segFirst, dcInfo->segLast );}



/***/
/* Select DCs in dTabAn1:
/***/

   validCountD  = 0;

   validCountD  = selectD(  dTabAn1->dObsList,
                            dTabAn1->dCount,
                            dcInfo->segFirst,
                            dcInfo->segLast,
                            dcInfo->excludeMode );


   validCountD  = selTypDC( dTabAn1->dObsList, dTabAn1->dCount, atomNameS1, atomNameS2 );

   FPR( stdout,
            "REMARK %d couplings selected (%s%s).\n",
            validCountD,
            dcInfo->excludeMode == DC_EXCLUDE     ? "exclusive" : "inclusive",
            dcInfo->excludeMode == DC_INCLUDE_MAX ? " max" : "" );



/***/
/* Compare DCs:
/***/

   validCountD  = adjustDC( palesInfo->opCode, palesInfo->c1,
                            palesInfo->c2, dTabAn1, dTabAn2 );

   FPR( stdout,
            "REMARK %d couplings combined (%s%s).\n",
            validCountD,
            dcInfo->excludeMode == DC_EXCLUDE     ? "exclusive" : "inclusive",
            dcInfo->excludeMode == DC_INCLUDE_MAX ? " max" : "" );



/***/
/* Output:
/***/

   (void) writeDCmod( outName, "w", dTabAn1, palesInfo, 
                 result );
 


/***/
/* Exit point:
/***/

shutdown:

    (void) deAlloc( "pales", dTabAn1,    sizeof( struct DTable ));
    (void) deAlloc( "pales", dTabAn2,    sizeof( struct DTable ));

 
    return (error); 
}
/* ========================================================================= */


/* ========================================================================= */
int adjustDC( int opCode, float c1, float c2,
              struct DTable *dTab1, struct DTable *dTab2 )
{
    int             i,j,ok,validCountD;
    float           d;
    struct  DInfo   *dList1, *dList2;



/***/
/* Initialize:
/***/

    i           = 0;
    j           = 0;
    ok          = 0;
    validCountD = 0;

    d           = 0.0;

    dList1       = (struct DInfo *)NULL;
    dList2       = (struct DInfo *)NULL;

    dList1     = dTab1->dObsList;
    dList2     = dTab1->dObsList;

    if (!dTab1)            return( 0 );
    if (!dTab2)            return( 0 );
    if (!dTab1->dCount)    return( 0 );
    if (!dTab2->dCount)    return( 0 );

    if (!dTab1->dObsList)  return( 0 );
    if (!dTab2->dObsList)  return( 0 );


    dList1     = dTab1->dObsList;
    dList2     = dTab2->dObsList;


    for( i = 0; i < dTab1->dCount; i++ )
      {
          ok = 0;

          if( dList1[i].d < 999.0 )
            {
              ok = 0;

              for( j = 0; j < dTab2->dCount; j++ )
                 {
                  if( (dList1[i].resIDI == dList2[j].resIDI) &&
                      (dList1[i].resIDJ == dList2[j].resIDJ) &&
                      !strcmp( dList1[i].resNameI, dList2[j].resNameI )  &&
                      !strcmp( dList1[i].resNameJ, dList2[j].resNameJ )  &&
                      !strcmp( dList1[i].atomNameI, dList2[j].atomNameI )  &&
                      !strcmp( dList1[i].atomNameI, dList2[j].atomNameI )  &&
                      !strcmp( dList1[i].atomNameJ, dList2[j].atomNameJ )  &&
                      !strcmp( dList1[i].chainNameI, dList2[j].chainNameI )  &&
                      !strcmp( dList1[i].chainNameJ, dList2[j].chainNameJ )  &&
                      !strcmp( dList1[i].segNameI, dList2[j].segNameI )  &&
                      !strcmp( dList1[i].segNameJ, dList2[j].segNameJ ) )
                    {
                       ok = 1;
                       d  = addDC( opCode, c1, c2, dList1[i].d, dList2[j].d );
                       continue;
                    }
                }

             if( ok )
               {
                  dTab1->dCalcList[i].d      = d;
                  dTab1->dCalcList[i].valid  = 1;
                  if( dList1[i].valid )
                    {
                       validCountD++;
                    }
              }
             else
               {
                  dList1[i].valid = 0;
                  dTab1->dCalcList[i].valid = 0;
               }
            }
       }



/***/
/* Exit point:
/***/

shutdown:

 
    return (validCountD); 
}
/* ========================================================================= */

/* ========================================================================= */
float addDC( int opCode, float c1, float c2,
                float d1, float d2 )
{
   float d;

   d = 0;


   switch( opCode )
      {
       case OP_ADD:
           d = (c1 * d1) + (c2 * d2);
           break;

       case OP_SUB:
           d = (c1 * d1) - (c2 * d2);
           break;

       case OP_MUL:
           d = (c1 * d1) * (c2 * d2);
           break;

       case OP_DIV:
           d = (c1 * d1) / (c2 * d2);
           break;

       case OP_MAX:
          if( d1 >= d2 )
            {
              d = c1 * d1;
            }
          else
            {
              d = c2 * d2;
            }
          break;

       case OP_MIN:
          if( d1 <= d2 )
            {
              d = c1 * d1;
            }
          else
            {
              d = c2 * d2;
            }
          break;
 
       default:
          break;
      }


    return( d );
}
/* ========================================================================= */

/* ========================================================================= */
int selTypDC( struct DInfo *dList, int dCount, char *atomNameI, 
              char *atomNameJ )
{
    int  i, validCount;


    validCount = 0;

    for( i = 0; i < dCount; i++ )
       {
        if ( !( !strcmp( atomNameI, "All" ) && !strcmp( atomNameJ, "All" ) ) &&
             !( (!strcmp( atomNameI, "All" ) && !strcmp( dList[i].atomNameJ, atomNameJ )) ||
                 (!strcmp( atomNameJ, "All" ) && !strcmp( dList[i].atomNameI, atomNameI )) ) &&
             !(!strcmp( dList[i].atomNameI, atomNameI ) && !strcmp( dList[i].atomNameJ, atomNameJ )) )
              dList[i].valid *= 0;
 
        if (dList[i].valid) validCount++;
       }

    return( validCount );
}
/* ========================================================================= */

