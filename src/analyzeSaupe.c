#include <stdio.h>
#include <string.h>
#include <math.h>


#define rnhn 	(1.042)
#define hCons	6.6260755e-4		/* as given in Bronstein */
#define hbCons	1.05457266e-4		/* as given in Bronstein */
#define u0	PI*4.0e-7
#define cx	((u0*hbCons)/(4.0*PI*PI))
#define gH	2.675198e+8
#define gN	-2.7116e+7


#define NO_LIST 4


#define PI 3.14159265

#define Real double


/* =================================================================== */
/* =================================================================== */
/* analyzeSaupe: analyze the Saupe matrix. 	       */
/* =================================================================== */
/* =================================================================== */
int analyzeSaupe(saupe, irrSaupe, genMagnitude,
        eigVec, eigVal, da, dr,
	alphalist, betalist, gammalist, verbose)

        int  verbose;
	Real *da,*dr, *genMagnitude;
	Real saupe[5], irrSaupe[5];
	Real eigVec[3][3], eigVal[3];
	Real alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
{
   Real da_in, dr_in, gMag;
   Real thetaGlist[NO_LIST], phiGlist[NO_LIST], omegaGlist[NO_LIST];
   Real psiList[8], thetaList[8], phiList[8];


   if( saupe2eulerMz( saupe, psiList, thetaList, phiList, 
	alphalist, betalist, gammalist, eigVec,
	eigVal, verbose) )
      printf("saupe2eulerMz: error occured\n");


   if(euler2polar(alphalist,betalist,gammalist,
	thetaGlist,phiGlist,omegaGlist))
      printf("euler2polar: error occured\n");


   if(pca_tensor(eigVal,&da_in,&dr_in))
      printf("pca_tensor: error occured\n");


/*
   if(print_analyzation( saupe, eigVec, eigVal, da_in, dr_in, 
	psiList, thetaList, phiList, 
	alphalist, betalist, gammalist,
	thetaGlist, phiGlist, omegaGlist))
      printf("print_analyzation: error occured\n");
*/


   irreducRepres( saupe, &gMag, irrSaupe );


   *genMagnitude = gMag;
   *da=da_in;
   *dr=dr_in;

   return 0;
}
/* =================================================================== */


/* =================================================================== */
int saupe2eulerMz( saupe, psiList, thetaList, phiList, 
	alphalist, betalist, gammalist, eigenvec,
	eigenvalue, verbose)

   int  verbose;
   Real alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
   Real saupe[5], psiList[8], thetaList[8], phiList[8];
   Real eigenvec[3][3], eigenvalue[3];
{
    int  i, j, n=0, nrot;
    Real sMat[4][4];
    Real eigVec[4][4], eigVal[4], eigVal_temp[4];


    for( i = 0; i < 8; i++ ) {
        psiList[i]   = 0.0;
        thetaList[i] = 0.0;
        phiList[i]   = 0.0;
    }

    for( i = 0; i < NO_LIST; i++ ) {
 	alphalist[i] = 0.0;
	betalist[i] = 0.0;
	gammalist[i] = 0.0;
    }

    sMat[1][1] = -0.5*(saupe[0] - saupe[1]);
    sMat[1][2] = saupe[2];
    sMat[1][3] = saupe[3];

    sMat[2][1] = saupe[2];
    sMat[2][2] = -0.5*(saupe[0] + saupe[1]);
    sMat[2][3] = saupe[4];

    sMat[3][1] = saupe[3];
    sMat[3][2] = saupe[4];
    sMat[3][3] = saupe[0];
 


   /*** determine eigensystem ***/
   jacobimz(sMat,eigVal,eigVec,&nrot);


   /*** sorting (absolute value) of eigenvalues and eigenvectors ***/
   if(sort_eigensystem(eigVal,eigVec))
	printf("sort_eigensystem: error occured\n");


   /*** determine Euler angles according to y-convention ***/
   /*** take into account ambiguity of eigenvectors *(+/-1) ***/
   if( verbose )
     {
       if(eigVec2euler(1,eigVec,alphalist,betalist,gammalist))
	  printf("eigVec2euler: error occured\n");
     }


   /*** determine Euler angles for rotation about three ***/
   /*** independent axes ***/
   /* not properly implemented */
   /* ambiguity of eigenvectors not taken into account */
   /* symmetry of model not taken into account */
/*
   for(i=1;i<=3;i++)
	for(j=1;j<=3;j++)	eigenvec[i-1][j-1] = eigVec[i][j];
   n = rot2euler(eigenvec,psiList,thetaList,phiList);
*/


   /*** return values */
   for(i=1;i<=3;i++){
	for(j=1;j<=3;j++)	eigenvec[i-1][j-1] = eigVec[i][j];
	eigenvalue[i-1] = eigVal[i];
   }


    return 0;
}
/* =================================================================== */




/* ===================== print alignment ========================= */
int print_analyzation( saupe, eigVec, eigVal, da, dr, 
	psiList, thetaList, phiList, 
	alphalist, betalist, gammalist,
	thetaGlist, phiGlist, omegaGlist)

	Real da, dr;
	Real saupe[5], eigVal[3], eigVec[3][3];
	Real psiList[8],thetaList[8],phiList[8];
	Real alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
	Real thetaGlist[NO_LIST], phiGlist[NO_LIST], omegaGlist[NO_LIST];
{
   int i,j;
   Real sMat[3][3];
   float daHz=0.0,R=0.0;

    sMat[0][0] = -0.5*(saupe[0] - saupe[1]);
    sMat[0][1] = saupe[2];
    sMat[0][2] = saupe[3];

    sMat[1][0] = saupe[2];
    sMat[1][1] = -0.5*(saupe[0] + saupe[1]);
    sMat[1][2] = saupe[4];

    sMat[2][0] = saupe[3];
    sMat[2][1] = saupe[4];
    sMat[2][2] = saupe[0];


    /*** print order matrix ***/
/*
    printf("Saupe matrix\n");
    for(i=0;i<3;i++){
	for(j=0;j<3;j++) printf("%e  ", sMat[i][j]);
	printf("\n");
    }
*/

   /*** print eigensystem ***/
/*
   printf("Eigenvalues\n");
   for(i=0;i<3;i++) printf("%e  ", eigVal[i]);
   printf("\n");
   printf("Eigenvectors\n");
    for(i=0;i<3;i++){
	for(j=0;j<3;j++) printf("%e  ", eigVec[i][j]);
	printf("\n");
    }
*/

    /*** print tensor characteristica ***/
    printf("\nMagnitude of alignment tensor\n");
    printf("==============================\n");
    daHz = da*(-1*gH*gN*cx/(rnhn*rnhn*rnhn));
    if(da!=0.0)	R = dr/da;
    printf("Da(NH) =  %.1f [Hz]     R = Dr/Da = %.3f\n\n",daHz,R);
    printf("\nOrientation of alignment tensor\n");
    printf("(Euler angles [deg.] for 4 equivalent orientations)\n");
    printf("==================================================\n");
    printf("alpha   :  ");
    for (i=0;i<NO_LIST;i++) printf("%7.2f  ",alphalist[i]);
    printf("\nbeta    :  ");
    for (i=0;i<NO_LIST;i++) printf("%7.2f  ",betalist[i]);
    printf("\ngamma   :  ");
    for (i=0;i<NO_LIST;i++) printf("%7.2f  ",gammalist[i]);
    printf("\n\n");

/*
    printf("Polar angles\n");
    printf("thetaGlist\n");
    for (i=0;i<NO_LIST;i++) printf("%6.2f  ",thetaGlist[i]);
    printf("\nphiGlist\n");
    for (i=0;i<NO_LIST;i++) printf("%6.2f  ",phiGlist[i]);
    printf("\nomegaGlist\n");
    for (i=0;i<NO_LIST;i++) printf("%6.2f  ",omegaGlist[i]);
    printf("\n\n");

    printf("\n\n");
    printf("\nEuler angles for rotation about three independent axes\n");
    printf("psiList\n");
    for (i=0;i<8;i++) printf("%.2f  ",psiList[i]);
    printf("\nthetaList\n");
    for (i=0;i<8;i++) printf("%.2f  ",thetaList[i]);
    printf("\nphiList\n");
    for (i=0;i<8;i++) printf("%.2f  ",phiList[i]);
    printf("\n\n");
*/


    return 0;
}
/* =============================================================== */


/* ========================================================================= */
int irreducRepres( saupe, genMagnitude, irrSaupe )
	double saupe[5], irrSaupe[5];
	double *genMagnitude;
{

    int i, error; 
    double fac, fac1, fac2, kappa;
    double szz, sxx_yy, sxy, sxz, syz, sxx, syy;


/***/
/* Initialize:
/***/

    error = 0;

    kappa  = 21523;
    fac    = (2.0/3.0) * (4.0/5.0*PI) * kappa*kappa;
    fac1   = sqrt(4.0/5.0*PI) * kappa;
    fac2   = sqrt(2.0/3.0 * 4.0/5.0*PI) * kappa;

    szz    = 0;
    szz    = 0;
    sxx_yy = 0;
    sxy    = 0;
    sxz    = 0;
    syz    = 0;
    sxx    = 0;
    syy    = 0;



/***/
/* Convert elements of Saupe matrix:
/***/

    szz    = saupe[0];
    sxx_yy = saupe[1];
    sxy    = saupe[2];
    sxz    = saupe[3];
    syz    = saupe[4];

    sxx    = -0.5*(szz - sxx_yy);
    syy    = -0.5*(szz + sxx_yy);



/***/
/* Calculate general magnitude of alignment tensor:
/***/

    *genMagnitude   = sqrt( fac *(sxx*sxx + syy*syy + szz*szz + 2*sxy*sxy + 2*sxz*sxz + 2*syz*syz ));

    irrSaupe[0]   = fac1 * szz;			/* a0 */
    irrSaupe[1]   = -1.0 * fac2 * sxz;  	/* a1r */
    irrSaupe[2]   = fac2 * syz;   		/* a1i */
    irrSaupe[3]   = fac2 * (sxx+szz/2.0);	/* a2r */
    irrSaupe[4]   = -1.0 * fac2 * sxy;		/* a2i */



/***/
/* Exit point:
/***/

shutdown:

    return (error); 
}
/* ========================================================================= */


/* ========================================================================= */
int irrCollinarity( saupe1, saupe2, dcos, angle )
	double saupe1[5], saupe2[5];
	float *dcos, *angle;
{

    int    i, error; 
    float  fac, fac1, fac2, kappa;
    float  szz1, sxx_yy1, sxy1, sxz1, syz1, sxx1, syy1;
    float  szz2, sxx_yy2, sxy2, sxz2, syz2, sxx2, syy2;
    float  gMag1, gMag2, coll, x;

/***/
/* Initialize:
/***/

    error = 0;

    kappa  = 21523;
    fac    = (2.0/3.0) * (4.0/5.0*PI) * kappa*kappa;
    fac1   = sqrt(4.0/5.0*PI) * kappa;
    fac2   = sqrt(2.0/3.0 * 4.0/5.0*PI) * kappa;

    szz1    = 0.0;
    szz1    = 0.0;
    sxx_yy1 = 0.0;
    sxy1    = 0.0;
    sxz1    = 0.0;
    syz1    = 0.0;
    sxx1    = 0.0;
    syy1    = 0.0;

    szz2    = 0.0;
    szz2    = 0.0;
    sxx_yy2 = 0.0;
    sxy2    = 0.0;
    sxz2    = 0.0;
    syz2    = 0.0;
    sxx2    = 0.0;
    syy2    = 0.0;

    gMag1   = 0.0;
    gMag2   = 0.0;
    coll    = 0.0;
    x       = 0.0;



/***/
/* Convert elements of Saupe matrix:
/***/

    szz1    = saupe1[0];
    sxx_yy1 = saupe1[1];
    sxy1    = saupe1[2];
    sxz1    = saupe1[3];
    syz1    = saupe1[4];

    sxx1    = -0.5*(szz1 - sxx_yy1);
    syy1    = -0.5*(szz1 + sxx_yy1);

    szz2    = saupe2[0];
    sxx_yy2 = saupe2[1];
    sxy2    = saupe2[2];
    sxz2    = saupe2[3];
    syz2    = saupe2[4];

    sxx2    = -0.5*(szz2 - sxx_yy2);
    syy2    = -0.5*(szz2 + sxx_yy2);


/***/
/* Calculate general magnitude of alignment tensor:
/***/

    coll    = fac * ( sxx1*sxx2 + syy1*syy2 + szz1*szz2 + 2*sxy1*sxy2 + 2*sxz1*sxz2 + 2*syz1*syz2 );

    gMag1   = sqrt( fac *(sxx1*sxx1 + syy1*syy1 + szz1*szz1 + 2*sxy1*sxy1 + 2*sxz1*sxz1 + 2*syz1*syz1 ));
    gMag2   = sqrt( fac *(sxx2*sxx2 + syy2*syy2 + szz2*szz2 + 2*sxy2*sxy2 + 2*sxz2*sxz2 + 2*syz2*syz2 ));

    if( gMag1*gMag2 != 0 )
      {
        *dcos   = coll / (gMag1*gMag2);
     }
    else
      {
	*dcos   = 0.0;
      }
 

    x = *dcos;

    if( x < 1 && x > -1 )
      {
	*angle =  90.0 - atan2(x,sqrt(1-x*x)) * (180.0/PI);
      }
    else
      {
	*angle =  0.0;
      }



/***/
/* Exit point:
/***/

shutdown:

    return (error); 
}
/* ========================================================================= */
