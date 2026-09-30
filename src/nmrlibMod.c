#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "pales.h"

#include "cmndargs.h"
#include "memory.h"
#include "prec.h"
#include "atof.h"

#include "hsymbol.h"
#include "pdb.h"
#include "pdbtab.h"
#include "dc.h"

#include "dynamo.h"

#define FPR             (void)fprintf


/* ========================================================================= */
int getDCTensorParams( int argc, char **argv, struct DTensor *tInfo )
{
    int i, loc1, loc2, loc, nextLoc, count;


    if (!tInfo) return( 0 );

/***/
/* Extract Tensor mode:
/***/

    (void) dTensorNull( tInfo );

    if (flagLoc( argc, argv, "-svd" ))     tInfo->tensorMode = DC_SVD;
    if (flagLoc( argc, argv, "-saupe" ))   tInfo->tensorMode = DC_SAUPE; 
    if (flagLoc( argc, argv, "-align" ))   tInfo->tensorMode = DC_SAUPE; 
    if (flagLoc( argc, argv, "-fixed" ))   tInfo->tensorMode = DC_FIXED_ALL;
    if (flagLoc( argc, argv, "-nofixed" )) tInfo->tensorMode = DC_FIXED_NONE;
    if (flagLoc( argc, argv, "-dadr" ))    tInfo->tensorMode = DC_FIXED_DADR;

    if (flagLoc( argc, argv, "-euler" ))   tInfo->eulerFlag  = 1;
    if (flagLoc( argc, argv, "-noeuler" )) tInfo->eulerFlag  = 0;

/***/
/* Extract DA and DR values if any:
/***/

    count = 0;

    count += fltArgD( argc, argv, "-da", &tInfo->da );
    count += fltArgD( argc, argv, "-dr", &tInfo->dr );

    if (count > 0 && tInfo->tensorMode == DC_SVD)
       {
        tInfo->tensorMode = DC_FIXED_DADR;
       }

    if (tInfo->tensorMode == DC_FIXED_DADR || tInfo->tensorMode == DC_FIXED_ALL)
       {
        tInfo->daMin = tInfo->da;
        tInfo->daMax = tInfo->da;
        tInfo->drMin = tInfo->dr;
        tInfo->drMax = tInfo->dr;
       }

    count += fltArgD( argc, argv, "-daMin",    &tInfo->daMin );
    count += fltArgD( argc, argv, "-daMax",    &tInfo->daMax );

    count += fltArgD( argc, argv, "-drMin",    &tInfo->drMin );
    count += fltArgD( argc, argv, "-drMax",    &tInfo->drMax );

/***/
/* Extract angular limits if any:
/***/

    count = 0;

    count += fltArgD( argc, argv, "-psi",       &tInfo->psi );
    count += fltArgD( argc, argv, "-theta",     &tInfo->theta );
    count += fltArgD( argc, argv, "-phi",       &tInfo->phi );

    if (tInfo->tensorMode == DC_FIXED_ALL)
       {
        tInfo->psiMin   = tInfo->psi;
        tInfo->psiMax   = tInfo->psi;

        tInfo->thetaMin = tInfo->theta;
        tInfo->thetaMax = tInfo->theta;

        tInfo->phiMin   = tInfo->phi;
        tInfo->phiMax   = tInfo->phi;
       }

    count += fltArgD( argc, argv, "-psiMin",    &tInfo->psiMin );
    count += fltArgD( argc, argv, "-psiMax",    &tInfo->psiMax );

    count += fltArgD( argc, argv, "-thetaMin",  &tInfo->thetaMin );
    count += fltArgD( argc, argv, "-thetaMax",  &tInfo->thetaMax );

    count += fltArgD( argc, argv, "-phiMin",    &tInfo->phiMin );
    count += fltArgD( argc, argv, "-phiMax",    &tInfo->phiMax );

/***/
/* Extract saupe values if any:
/***/

    loc  = 0;
    loc1 = flagLoc( argc, argv, "-align" );
    loc2 = flagLoc( argc, argv, "-saupe" );

    if (loc1) loc = loc1;
    if (loc2) loc = loc2;

    if (loc) 
       {
        nextLoc = nextFlag( argc, argv, loc );

        if (nextLoc - loc != 6)
           {
            FPR( stderr, "DYNAMO Error: -saupe option requires 5 values.\n" );
            return( 1 );
           }

        for( i = 0; i < 5; i++ )
           {
            tInfo->userSaupe[i] = ATOF( getNthArg( argc, argv, loc+i+1 ));
           }
       }

    return( 0 );
}
/* ========================================================================= */

/* ========================================================================= */
int writeDCmod( outDName, access, dTab, palesInfo, remark )

   struct DTable  *dTab;
   struct PalesInfo *palesInfo;
   char   *outDName, *access, *remark;
{
    struct  DInfo   *dObsPtr, *dCalcPtr;
    struct  DTensor *dTensor;

    float   d, rSumD2,pSumD,qD;
    float   coorR,offset,errOffset;
    float   slope,errSlope,aveSlope,errAveSlope;

    double  saupe[5], irrSaupe[5];
    double  eigVal[3],eigVec[3][3];
    double  gMag,da,dr,R;
    static  double alphaList[NO_LIST], betaList[NO_LIST], gammaList[NO_LIST];

    int     i, error, nD;

    FILE    *outUnit;


/***/
/* Initialize:
/***/

    i           = 0;
    error       = 0;
    rSumD2      = 0.0;
    pSumD       = 0.0;
    qD          = 0.0;
    coorR       = 0.0;
    offset      = 0.0;
    errOffset   = 0.0;
    slope       = 0.0;
    errSlope    = 0.0;
    aveSlope    = 0.0;
    errAveSlope = 0.0;
    nD          = 0;

    dTensor   = &dTab->dTensor;


/***/
/* Open the output file:
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

    /* Statistics and nD are needed by both verbose and compact output. */
    (void) getStatDCmod( dTab,
        &rSumD2, &pSumD, &qD, &coorR,
        &offset, &errOffset, &slope, &errSlope,
        &aveSlope, &errAveSlope, &nD );


 if (palesInfo->verbFlag)
   {


/***/
/* Header:
/***/

    FPR( outUnit, "REMARK Molecular Alignment Simulation.\n\n" );

    if (remark) FPR( outUnit, "%s", remark );


    FPR( outUnit,
          "REMARK Simulation parameters.\n" );


/***/
/* Report pales info:
/***/

    switch( palesInfo->palesMode )
      {
       case PALES_STERIC:
           FPR( outUnit, "\nDATA PALES_MODE STERIC\n\n" );
  
           if( palesInfo->lcType == 1 )
             {
               FPR( outUnit, "DATA PALES LC_TYPE             wall\n");        
             }
           else if( palesInfo->lcType == 2 )
             {
              FPR( outUnit, "DATA PALES LC_TYPE              rod\n" );        
             }
 
            FPR( outUnit, "DATA PALES LC_CONCENTRATION    %.3f\n", palesInfo->concentration);        
            FPR( outUnit, "DATA PALES ORIENT_SPHERE       %d\n", palesInfo->densityNSC);        
            FPR( outUnit, "DATA PALES ORIENT_PSI          %d\n", palesInfo->digPsiNSC);        
            FPR( outUnit, "DATA PALES GRID_SPACING       % .3f\n", palesInfo->gridSpacing);        
            FPR( outUnit, "DATA PALES MODEL_RADIUS       % .3f\n", palesInfo->radiusModel);        
            FPR( outUnit, "DATA PALES LC_ORDER           % .3f\n", palesInfo->lcOrder);        
            FPR( outUnit, "DATA PALES ATOM_RADIUS        % .3f\n", palesInfo->radiusAtom);        
            FPR( outUnit, "DATA PALES SEL_SIMPLE_FLAG     %d\n", palesInfo->Hflag);        
            FPR( outUnit, "DATA PALES SURF_FLAG           %d\n", palesInfo->surfFlag);        

           break;

       case PALES_STERIC_FREE:
           FPR( outUnit, "\nDATA PALES_MODE STERIC_FREE\n\n" );

           if( palesInfo->lcType == 1 )
             {
               FPR( outUnit, "DATA PALES LC_TYPE             wall\n");        
             }
           else if( palesInfo->lcType == 2 )
             {
              FPR( outUnit, "DATA PALES LC_TYPE              rod\n" );        
             }
 
            FPR( outUnit, "DATA PALES LC_CONCENTRATION    %.3f\n", palesInfo->concentration);        
            FPR( outUnit, "DATA PALES ORIENT_SPHERE       %d\n", palesInfo->densityNSC);        
            FPR( outUnit, "DATA PALES ORIENT_PSI          %d\n", palesInfo->digPsiNSC);        
            FPR( outUnit, "DATA PALES GRID_SPACING       % .3f\n", palesInfo->gridSpacing);        
            FPR( outUnit, "DATA PALES MODEL_RADIUS       % .3f\n", palesInfo->radiusModel);        
            FPR( outUnit, "DATA PALES LC_ORDER           % .3f\n", palesInfo->lcOrder);        
            FPR( outUnit, "DATA PALES ATOM_RADIUS        % .3f\n", palesInfo->radiusAtom);        
            FPR( outUnit, "DATA PALES SEL_SIMPLE_FLAG     %d\n", palesInfo->Hflag);        
            FPR( outUnit, "DATA PALES SURF_FLAG           %d\n", palesInfo->surfFlag);        

           break;

       case PALES_STERIC_GRID:
            FPR( outUnit, "\nDATA PALES_MODE STERIC_GRID\n\n" );

           if( palesInfo->lcType == 1 )
             {
               FPR( outUnit, "DATA PALES LC_TYPE             wall\n");        
             }
           else if( palesInfo->lcType == 2 )
             {
              FPR( outUnit, "DATA PALES LC_TYPE              rod\n" );        
             }
 
            FPR( outUnit, "DATA PALES LC_CONCENTRATION    %.3f\n", palesInfo->concentration);        
            FPR( outUnit, "DATA PALES ORIENT_SPHERE       %d\n", palesInfo->densityNSC);        
            FPR( outUnit, "DATA PALES ORIENT_PSI          %d\n", palesInfo->digPsiNSC);        
            FPR( outUnit, "DATA PALES GRID_SPACING       % .3f\n", palesInfo->gridSpacing);        
            FPR( outUnit, "DATA PALES MODEL_RADIUS       % .3f\n", palesInfo->radiusModel);        
            FPR( outUnit, "DATA PALES LC_ORDER           % .3f\n", palesInfo->lcOrder);        
            FPR( outUnit, "DATA PALES ATOM_RADIUS        % .3f\n", palesInfo->radiusAtom);        
            FPR( outUnit, "DATA PALES SEL_SIMPLE_FLAG     %d\n", palesInfo->Hflag);        
            FPR( outUnit, "DATA PALES SURF_FLAG           %d\n", palesInfo->surfFlag);        

           break;

       case PALES_ELECTROSTATIC:
            FPR( outUnit, "\nDATA PALES_MODE ELECTROSTATIC\n\n" );

           if( palesInfo->lcType == 1 )
             {
               FPR( outUnit, "DATA PALES LC_TYPE             wall\n");        
             }
           else if( palesInfo->lcType == 2 )
             {
              FPR( outUnit, "DATA PALES LC_TYPE              rod\n" );        
             }
 
            FPR( outUnit, "DATA PALES LC_CONCENTRATION    %.3f\n", palesInfo->concentration);        
            FPR( outUnit, "DATA PALES ORIENT_SPHERE       %d\n", palesInfo->densityNSC);        
            FPR( outUnit, "DATA PALES ORIENT_PSI          %d\n", palesInfo->digPsiNSC);        
            FPR( outUnit, "DATA PALES GRID_SPACING       % .3f\n", palesInfo->gridSpacing);        
            FPR( outUnit, "DATA PALES MODEL_RADIUS       % .3f\n", palesInfo->radiusModel);        
            FPR( outUnit, "DATA PALES LC_ORDER           % .3f\n", palesInfo->lcOrder);        
            FPR( outUnit, "DATA PALES ATOM_RADIUS        % .3f\n", palesInfo->radiusAtom);        
            FPR( outUnit, "DATA PALES SEL_SIMPLE_FLAG     %d\n", palesInfo->Hflag);        
            FPR( outUnit, "DATA PALES SURF_FLAG           %d\n", palesInfo->surfFlag);        

            FPR( outUnit, "DATA PALES EL PH              % .3f\n", palesInfo->elInfo.pH);
            FPR( outUnit, "DATA PALES EL NACL            % .3f\n", palesInfo->elInfo.concNaCl);
            FPR( outUnit, "DATA PALES EL SURF_CH         % .3f\n", palesInfo->elInfo.chDensSurf);

            if( palesInfo->elInfo.pkaFlag )
             {
               FPR( outUnit, "DATA PALES EL PKA_DEFAULT         on\n" );
             }
           else
             {
               FPR( outUnit, "DATA PALES EL PKA_DEFAULT         off\n" );
             }
     
            if( palesInfo->elInfo.chTerminiFlag )
             {
               FPR( outUnit, "DATA PALES EL CH_TERM             on\n" );
             }
           else
             {
               FPR( outUnit, "DATA PALES EL CH_TERM             off\n" );
             }

            FPR( outUnit, "DATA PALES EL CUTOFF         % .3f\n", palesInfo->elInfo.cutoffPot);

           break;

       case PALES_ELECTROSTATIC_FREE:
            FPR( outUnit, "\nDATA PALES_MODE ELECTROSTATIC_FREE\n\n" );

           if( palesInfo->lcType == 1 )
             {
               FPR( outUnit, "DATA PALES LC_TYPE             wall\n");        
             }
           else if( palesInfo->lcType == 2 )
             {
              FPR( outUnit, "DATA PALES LC_TYPE              rod\n" );        
             }
 
           FPR( outUnit, "DATA PALES LC_CONCENTRATION    %.3f\n", palesInfo->concentration);      
           FPR( outUnit, "DATA PALES ORIENT_SPHERE       %d\n", palesInfo->densityNSC);        
           FPR( outUnit, "DATA PALES ORIENT_PSI          %d\n", palesInfo->digPsiNSC);        
           FPR( outUnit, "DATA PALES GRID_SPACING       % .3f\n", palesInfo->gridSpacing);        
           FPR( outUnit, "DATA PALES MODEL_RADIUS       % .3f\n", palesInfo->radiusModel);        
           FPR( outUnit, "DATA PALES LC_ORDER           % .3f\n", palesInfo->lcOrder);        
           FPR( outUnit, "DATA PALES ATOM_RADIUS        % .3f\n", palesInfo->radiusAtom);        
           FPR( outUnit, "DATA PALES SEL_SIMPLE_FLAG     %d\n", palesInfo->Hflag);        
           FPR( outUnit, "DATA PALES SURF_FLAG           %d\n", palesInfo->surfFlag);        

           break;
 
       case PALES_ELECTROSTATIC_GRID:
           FPR( outUnit, "\nDATA PALES_MODE ELECTROSTATIC_GRID\n\n" );

           if( palesInfo->lcType == 1 )
             {
               FPR( outUnit, "DATA PALES LC_TYPE             wall\n");        
             }
           else if( palesInfo->lcType == 2 )
             {
              FPR( outUnit, "DATA PALES LC_TYPE              rod\n" );        
             }
 
           FPR( outUnit, "DATA PALES LC_CONCENTRATION    %.3f\n", palesInfo->concentration);        
           FPR( outUnit, "DATA PALES ORIENT_SPHERE       %d\n", palesInfo->densityNSC);        

           FPR( outUnit, "DATA PALES ORIENT_PSI          %d\n", palesInfo->digPsiNSC);        
           FPR( outUnit, "DATA PALES GRID_SPACING       % .3f\n", palesInfo->gridSpacing);        
           FPR( outUnit, "DATA PALES MODEL_RADIUS       % .3f\n", palesInfo->radiusModel);        
           FPR( outUnit, "DATA PALES LC_ORDER           % .3f\n", palesInfo->lcOrder);        
           FPR( outUnit, "DATA PALES ATOM_RADIUS        % .3f\n", palesInfo->radiusAtom);        
           FPR( outUnit, "DATA PALES SEL_SIMPLE_FLAG     %d\n", palesInfo->Hflag);        
           FPR( outUnit, "DATA PALES SURF_FLAG           %d\n", palesInfo->surfFlag);        

           break;
 
       case PALES_DC:
            FPR( outUnit, "\nDATA PALES_MODE DC\n\n" );

           switch( dTensor->tensorMode )
             {
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
                   dTensor->psiMin, dTensor->psiMax );

                 FPR( outUnit,
                   "DATA THETA_RANGE % 8.2f % 8.2f\n",
                   dTensor->thetaMin, dTensor->thetaMax );

                 FPR( outUnit,
                   "DATA PHI_RANGE   % 8.2f % 8.2f\n",
                   dTensor->phiMin, dTensor->phiMax );
                 break;

               case DC_FIXED_NONE:
                 FPR( outUnit, "DATA TENSOR_MODE NOFIXED" );
                 FPR( outUnit, " (Variable Da and Dr, Variable Rotations)\n\n" );

                 FPR( outUnit,
                   "DATA DA_RANGE    % .4e % .4e\n",
                   dTensor->daMin, dTensor->daMax );

                 FPR( outUnit,
                   "DATA DR_RANGE    % .4e % .4e\n",
                   dTensor->drMin, dTensor->drMax );

                 FPR( outUnit,
                   "DATA PSI_RANGE   % 8.2f % 8.2f\n",
                   dTensor->psiMin, dTensor->psiMax );

                 FPR( outUnit,
                   "DATA THETA_RANGE % 8.2f % 8.2f\n",
                   dTensor->thetaMin, dTensor->thetaMax );

                 FPR( outUnit,
                   "DATA PHI_RANGE   % 8.2f % 8.2f\n",
                   dTensor->phiMin, dTensor->phiMax );
                 break;

               default:
                 break;
           }

           break;
 
       case PALES_CS:
           FPR( outUnit, "\nDATA PALES_MODE CS\n\n" );
           break;
 
       case PALES_HIST:
           FPR( outUnit, "\nDATA PALES_MODE HIST\n\n" );
           break;

       case PALES_ASSIGN:
           FPR( outUnit, "\nDATA PALES_MODE HIST\n\n" );
           break;

       case PALES_DC_CP:

           FPR( outUnit, "\nDATA PALES_MODE DC_CP\n\n" );

           switch( palesInfo->opCode )
             {
              case OP_ADD:
                FPR( outUnit, "DATA PALES ACTION    add\n" );
                 break;

              case OP_SUB:
                 FPR( outUnit, "DATA PALES ACTION    sub\n" );
                 break;

             case OP_MUL:
                 FPR( outUnit, "DATA PALES ACTION    mul\n" );
                 break;

            case OP_DIV:
                 FPR( outUnit, "DATA PALES ACTION    div\n" );
                 break;

            case OP_MAX:
                 FPR( outUnit, "DATA PALES ACTION    max\n" );
                 break;
            default:
                 break;
            }
 
           FPR( outUnit, "DATA PALES C1     % .3f\n", palesInfo->c1 );
           FPR( outUnit, "DATA PALES C2     % .3f\n", palesInfo->c2 );

           break;

       case PALES_SAUPE_CP:
           FPR( outUnit, "\nDATA PALES_MODE SAUPE_CP\n\n" );

           switch( palesInfo->opCode )
             {
              case OP_ADD:
                FPR( outUnit, "DATA PALES ACTION    add\n" );
                 break;

              case OP_SUB:
                 FPR( outUnit, "DATA PALES ACTION    sub\n" );
                 break;

             case OP_MUL:
                 FPR( outUnit, "DATA PALES ACTION    mul\n" );
                 break;

            case OP_DIV:
                 FPR( outUnit, "DATA PALES ACTION    div\n" );
                 break;

            case OP_MAX:
                 FPR( outUnit, "DATA PALES ACTION    max\n" );
                 break;
            default:
                 break;
            }
 
           FPR( outUnit, "DATA PALES C1     % .3f\n", palesInfo->c1 );
           FPR( outUnit, "DATA PALES C2     % .3f\n", palesInfo->c2 );

           break;

       default:
           break;
      }


    if( palesInfo->includeRim != 0 )
      {
        FPR( outUnit, "DATA PALES INCLUDE_RIM      on\n" );
        FPR( outUnit, "DATA PALES MOLAR_RATIO      % .3f\n", palesInfo->molarRatio );
        FPR( outUnit, "DATA PALES MODEL_LENGTH     % .3f\n", palesInfo->lengthModel );
        FPR( outUnit, "DATA PALES S_DISC_RIM       % .3f\n", palesInfo->orderDiscRim );
      }


    FPR( outUnit, "\n\n" );



/***/
/* Report saupe matrix:
/***/
    
    FPR( outUnit,
          "REMARK Order matrix.\n\n" );

    FPR( outUnit,
           "DATA SAUPE  % .4e % .4e % .4e % .4e % .4e\n\n",
           dTensor->saupe[0], 
           dTensor->saupe[1], 
           dTensor->saupe[2], 
           dTensor->saupe[3], 
           dTensor->saupe[4] );

 
    saupe[0] = dTensor->saupe[0]; 
    saupe[1] = dTensor->saupe[1]; 
    saupe[2] = dTensor->saupe[2]; 
    saupe[3] = dTensor->saupe[3]; 
    saupe[4] = dTensor->saupe[4];



/***/
/* Analyze order matrix:
/***/

   if( saupe[0] != 0 &&
       saupe[1] != 0 &&
       saupe[2] != 0 &&
       saupe[3] != 0 &&
       saupe[4] != 0 )
     {
       if( analyzeSaupe(saupe, irrSaupe, &gMag, eigVec, eigVal, &da, &dr,
           alphaList, betaList, gammaList, palesInfo->eulerFlag) )
         {
           FPR( stderr, "PALES: error occured in analyzeSaupe.\n");
           error = 1;
           goto shutdown;
         }


 
/***/
/* Write eigemvalues & eigenvectors:
/***/

       FPR( outUnit,
             "DATA IRREDUCIBLE REPRESENTATION (A0,A1R,A1I,A2R,A2I)" );

       FPR( outUnit,
             "   % .4e % .4e % .4e % .4e % .4e\n",
             irrSaupe[0], irrSaupe[1], irrSaupe[2], irrSaupe[3], irrSaupe[4] );


       FPR( outUnit,
             "DATA IRREDUCIBLE GENERAL_MAGNITUDE" );

       FPR( outUnit,
             "  % .4e\n",  gMag );


       if ( palesInfo->eulerFlag )
          {

             FPR( outUnit,
                    "\n\nREMARK Eigensystem & Euler angles for clockwise rotation about z, y', z''.\n\n" );


             FPR( outUnit,
                    "DATA EIGENVALUES (Axx,Ayy,Azz)" );

             FPR( outUnit,
                    "   % .4e % .4e % .4e\n",
                    eigVal[2], eigVal[1], eigVal[0] );


             FPR( outUnit,
                    "DATA EIGENVECTORS\n" );

             FPR( outUnit,
                    "DATA EIGENVECTORS XAXIS % .4e % .4e % .4e\n",
                    eigVec[0][2],
                    eigVec[1][2],
                    eigVec[2][2] );

             FPR( outUnit,
                    "DATA EIGENVECTORS YAXIS % .4e % .4e % .4e\n",
                    eigVec[0][1],
                    eigVec[1][1],
                    eigVec[2][1] );

             FPR( outUnit,
                    "DATA EIGENVECTORS ZAXIS % .4e % .4e % .4e\n",
                    eigVec[0][0],
                    eigVec[1][0],
                    eigVec[2][0] );



/***/
/* Write Euler angles:
/***/

             FPR( outUnit, "\n" );

             FPR( outUnit,
                    "DATA Q_EULER_SOLUTIONS  %7s  %7s  %7s\n",
                    "ALPHA", "BETA", "GAMMA" );


             for( i = 0; i < NO_LIST; i++ )
               {
                  FPR( outUnit,
                  "DATA Q_EULER_ANGLES %2d  %7.2f  %7.2f  %7.2f\n",
                  i+1,
                  alphaList[i],
                  betaList[i],
                  gammaList[i] );
               }

             FPR( outUnit,
                  "\n\nREMARK Euler angles (psi/theta/phi) for rotation about x, y, z.\n\n" );

             (void) reportEuler( outUnit, dTensor );


          } /* end of if palesInfo->verbFlag */


     } /* end of euler output */


/***/
/* Calculate statistics of DCs:
/***/

/***/
/* Report statisctics of DCs:
/***/

    FPR( outUnit,
          "\nREMARK Dipolar couplings.\n\n" );

    FPR( outUnit, "DATA N                    %d\n",       nD );


   } /* end of if palesInfo->verbFlag */
 else
   {
    FPR( outUnit, "DATA SEQUENCE " );

/*
    for( i = 0; i < palesInfo->resCount; i++ )
      {
        j++;

        FPR( outUnit, "%s", palesInfo->res[i] );

        if( j == 10)
          { 
            FPR( outUnit, " " );
            j = 0; 
            k++; 
          }    
        if( k == 4 )   
          { 
            k = 0;
            FPR( outUnit, "\nDATA SEQUENCE " );
          }    
      }
*/

    FPR( outUnit, "\n" );
   }  


 if( nD != 0 )
   {

    if( palesInfo->verbFlag )
      {
        FPR( outUnit, "DATA RMS                  %.3f\n",     rSumD2 );
        FPR( outUnit, "DATA Chi2                 %.3f\n",     pSumD );
        FPR( outUnit, "DATA CORR R               %.3f\n",     coorR );
        FPR( outUnit, "DATA CORNILESCU Q         %.3f\n",     qD );
        FPR( outUnit, "DATA REGRESSION OFFSET    %.3f +/- %.3f [Hz]\n",
               offset, errOffset );
        FPR( outUnit, "DATA REGRESSION SLOPE     %.3f +/- %.3f [Hz]\n",
           slope, errSlope );
        FPR( outUnit, "DATA REGRESSION BAX SLOPE %.3f +/- %.3f [Hz]\n",
               aveSlope, errAveSlope );
      }
 

    FPR( outUnit, "\n" );



/***/
/* Report of simulated couplings:
/***/

    FPR( outUnit, "VARS   " );

    PrintAtomIVars( outUnit, dTab );
    PrintAtomJVars( outUnit, dTab );

    if (palesInfo->verbFlag)
      FPR( outUnit, " DI D_OBS D D_DIFF DD W\n" );
    else
      FPR( outUnit, " D DD W\n" );

    FPR( outUnit, "FORMAT " );
   
    PrintAtomIFmt( outUnit, dTab );
    PrintAtomJFmt( outUnit, dTab );

    if (palesInfo->verbFlag)
      FPR( outUnit, " %%9.2f %%9.3f %%9.3f %%9.3f %%.2f %%.2f\n\n" );
    else
      FPR( outUnit, " %%9.3f %%9.3f %%.2f\n\n" );

    dObsPtr  = dTab->dObsList;
    dCalcPtr = dTab->dCalcList;

 
    for( i = 0; i < dTab->dCount; i++ )
      {
        if (dCalcPtr->valid && dObsPtr->valid)
          {
            if (dObsPtr->d < 999.0 && dCalcPtr->d < 999.0)
               d = dObsPtr->d - dCalcPtr->d;
            else
               d = DC_UVAL;

            PrintAtomIVals( outUnit, dObsPtr, dTab );
            PrintAtomJVals( outUnit, dObsPtr, dTab );

            if( palesInfo->verbFlag )
              { 
                FPR( outUnit,
                  " % 9.2f % 9.4f % 9.4f % 9.4f %7.4f %.2f\n",
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
                  " % 9.4f %7.4f %.2f\n",
                  dCalcPtr->d,
                  dObsPtr->sigma,
                  dObsPtr->w );
              }
          }
 
       dObsPtr++; dCalcPtr++;
      }

   } /* end of if palesInfo->verbFlag && nD != 0 */



/***/
/* Exit point:
/***/

shutdown:

   if (outUnit != stdout) (void) fclose( outUnit );

   return( error );
}
/* ===================================================================== */

/* ========== compare measured with predicted dipolar couplings ======== */
int getStatDCmod(struct DTable *dTab, 
        float *rSumD2, float *pSumD, float *qD, float *coorR,
        float *offset,  float *errOffset, float *slope, float *errSlope,
        float *aveSlope, float *errAveSlope,
        int *nD )
{
   int    i;
   float  prob_stat, z_stat;
   float  *dcNorm,*dcMatrixNorm, *sigma; 
   float  *dcNorm2,*dcMatrixNorm2; 
   float  stdDevD, offNull,errOffNull,chi2,q;
   float  slope1, slope2, sigb1, sigb2, siga;
   float  rSumD,sumD2,d;

   int             error, dcCount;
   struct DInfo    *dList;
 


/***/
/* Initialize:
/***/

    error       = 0;
    i           = 0;

    dcCount     = 0;
    d           = 0;
    stdDevD     = 0;
    sumD2       = 0.0;
    rSumD       = 0.0;
    offNull     = 0;
    errOffNull  = 0;
    chi2        = 0;
    q           = 0;


    siga        = 0;
    slope1      = 0;
    slope2      = 0;
    sigb1       = 0;
    sigb2       = 0;

    dList       = (struct DInfo *)NULL;

    dcNorm        = (float *)NULL;
    dcMatrixNorm  = (float *)NULL;
    dcNorm2       = (float *)NULL;
    dcMatrixNorm2 = (float *)NULL;
    sigma         = (float *)NULL;

    dcCount   = dTab->dCount;
    dList     = dTab->dObsList;

/***/
/* Allocate and extract coord lists, etc:
/***/


    dcNorm        = fltAlloc( "statDC", dcCount+1 );
    dcMatrixNorm  = fltAlloc( "statDC", dcCount+1 );
    sigma         = fltAlloc( "statDC", dcCount+1 );
    dcNorm2       = fltAlloc( "statDC", dcCount+5001 );
    dcMatrixNorm2 = fltAlloc( "statDC", dcCount+5001 );

    if (!dcNorm || !dcMatrixNorm || !sigma || !dcNorm2 || !dcMatrixNorm2 )
       {
        FPR( stderr, "Error allocating temporary buffers.\n" );
        error = 1;
        goto shutdown;
       }


/***/
/* Do the statistics
/***/

    for( i = 0; i < dcCount; i++ )
       {
          if (!dList[i].valid) continue;

          if (dTab->dCalcList[i].valid && 
              dTab->dObsList[i].d < 999.0 && 
              dTab->dCalcList[i].d < 999.0)
           {
             (*nD)++;

              dcNorm[*nD]        = dTab->dObsList[i].d  * dTab->dCalcList[i].w;
              dcMatrixNorm[*nD]  = dTab->dCalcList[i].d * dTab->dCalcList[i].w;
              sigma[*nD]         = dTab->dObsList[i].sigma;

              d        = dcNorm[*nD] - dcMatrixNorm[*nD];

              sumD2   += dcNorm[*nD]*dcNorm[*nD];
              rSumD   += d;
              *rSumD2 += d*d;

              d       /= dTab->dObsList[i].sigma == 0.0 ? 1.0 : dTab->dObsList[i].sigma;
              *pSumD  += d*d;
           }
       }


     stdDevD = *nD ? sqrt( (double)((*rSumD2 - rSumD*rSumD / (*nD)) / (*nD)) ): 0.0;

     *qD      = sumD2 ? sqrt( (double)((*rSumD2)/(sumD2)) ) : 0.0;

     *rSumD2  = *nD ? sqrt( (double)((*rSumD2)/(*nD)) ) : 0.0;



     if( *nD != 0 )
       {
          /* linear correlation */
          pearsn(dcNorm,dcMatrixNorm,*nD,coorR,
		&prob_stat,&z_stat);

          /* linear regression */
          fit(dcNorm,dcMatrixNorm,*nD,sigma,0,
		offset,slope,errOffset,errSlope,&chi2,&q);
       }



/***/
/* Calculate average slope with offset = 0:
/***/

     if( *nD != 0 )
       {
          *aveSlope = 0;

          for( i = 1; i <= *nD; i++ )
             {
               dcNorm2[i]       = dcNorm[i];
               dcMatrixNorm2[i] = dcMatrixNorm[i];
             }

          for( i = (*nD)+1; i <= (*nD)+5000; i++ )
             {
               dcNorm2[i]       = 0;
               dcMatrixNorm2[i] = 0;
             }


          fit(dcNorm2,dcMatrixNorm2,*nD,sigma,0,
		&offNull,&slope1,&siga,&sigb1,&chi2,&q);

          fit(dcMatrixNorm2,dcNorm2,*nD,sigma,0,
		&offNull,&slope2,&siga,&sigb2,&chi2,&q);


          *aveSlope    = (slope2 != 0) ? 0.5*(slope1 + 1.0/slope2) : 0.0;

          *errAveSlope = (slope2 != 0) ? 0.5*(*aveSlope)* sqrt( sigb1*sigb1 + (sigb2/(slope2*slope2))*(sigb2/(slope2*slope2)) ) : 0.0;
       }



/***/
/* Exit point:
/***/

shutdown:

    (void) deAlloc( "statDC", dcNorm,        sizeof(float)*(dcCount+1) );
    (void) deAlloc( "statDC", dcMatrixNorm,  sizeof(float)*(dcCount+1) );
    (void) deAlloc( "statDC", sigma,         sizeof(float)*(dcCount+1) );
    (void) deAlloc( "statDC", dcNorm2,       sizeof(float)*(dcCount+5001) );
    (void) deAlloc( "statDC", dcMatrixNorm2, sizeof(float)*(dcCount+5001) );

    return( error );
}
/* ========================================================================= */

/* ========================================================================= */
int writeDcMz( outDName, access, dTab, remark, verbose )

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
int writeDcMzOld( outDName, access, dTab, remark, verbose )

   struct DTable  *dTab;
   char   *outDName, *remark, *access;
   int    verbose;
{
#define varStr1  "RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J"
#define fmtStr1  "%5d     %6s       %6s        %5d     %6s       %6s       "

#define varStr2  "DI D DD W"
#define fmtStr2  "%9.2f %9.3f %9.3f %.2f"

#define varStr2V "DI    D_OBS           D     D_DIFF DD   W"
#define fmtStr2V "%9.2f %9.3f           %9.3f %9.3f  %.2f %.2f"

    struct DInfo *dObsList, *dCalcList;
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

    if (verbose)
       {
        FPR( outUnit, "VARS   %s %s\n\n", varStr1, varStr2V );
        FPR( outUnit, "FORMAT %s %s\n\n", fmtStr1, fmtStr2V );
       }
    else
       {
        FPR( outUnit, "VARS   %s %s\n",   varStr1, varStr2 );
        FPR( outUnit, "FORMAT %s %s\n\n", fmtStr1, fmtStr2 );
       }

    for( i = 0; i < dTab->dCount; i++ )
       {
        if (dCalcList[i].valid)
           {
            if (dObsList[i].d < 999.0 && dCalcList[i].d < 999.0)
               d = dObsList[i].d - dCalcList[i].d;
            else
               d = DC_UVAL;

            if (verbose)
               {
                FPR( outUnit,
                     "%3d %4s %4s %3d %4s %4s ",
                     dObsList[i].resIDI,
                     dObsList[i].resNameI,
                     dObsList[i].atomNameI,
                     dObsList[i].resIDJ,
                     dObsList[i].resNameJ,
                     dObsList[i].atomNameJ );

                FPR( outUnit,
                     "% 9.2f % 9.4f % 9.4f % 9.4f %7.4f %.2f\n",
                     dCalcList[i].di,
                     dObsList[i].d,
                     dCalcList[i].d,
                     d,
                     dObsList[i].sigma,
                     dObsList[i].w );
               }
            else
               {
                FPR( outUnit,
                     "%4d %4s %4s %4d %4s %4s % 9.2f % 9.4f %7.4f %.2f\n",
                     dObsList[i].resIDI,
                     dObsList[i].resNameI,
                     dObsList[i].atomNameI,
                     dObsList[i].resIDJ,
                     dObsList[i].resNameJ,
                     dObsList[i].atomNameJ,
                     dCalcList[i].di,
                     dCalcList[i].d,
                     dObsList[i].sigma,
                     dObsList[i].w );
               }
           }
       }

    if (outUnit != stdout) (void) fclose( outUnit );

    return( 0 );
}
/* ========================================================================= */

/* ========================================================================= */
int getEulerDC (float eigenVec[3][3], struct DTensor *dTensor)
{
 
    dTensor->eulerCount = rot2euler( eigenVec, 
                                     dTensor->psiList, 
                                     dTensor->thetaList, 
                                     dTensor->phiList );

    (void) euler2saupe( dTensor->saupe,
                        dTensor->psiList[0], 
                        dTensor->thetaList[0], 
                        dTensor->phiList[0],
                        dTensor->da,
                        dTensor->dr );

    return( dTensor->eulerCount );
}
/* ========================================================================= */
