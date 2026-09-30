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

#define RSWAP( RA, RB )  rtemp = RA; RA = RB; RB = rtemp

#define RAD(x) (x*PI/180.)
#define DEG(x) (x*180./PI)



/* =================================================================== */
/* =================================================================== */
/* analyzeMatrix3D: analyze the Saupe matrix. 	       */
/* =================================================================== */
/* =================================================================== */
int analyzeMatrix3D(matrix, eigVec, eigVal, da, dr,
	alphalist, betalist, gammalist)

	Real *da,*dr;
	Real matrix[9];
	Real eigVec[3][3], eigVal[3];
	Real alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
{
   Real da_in,dr_in;
   Real thetaGlist[NO_LIST], phiGlist[NO_LIST], omegaGlist[NO_LIST];
   Real psiList[8], thetaList[8], phiList[8];


   if(matrix2euler( matrix, psiList, thetaList, phiList, 
	alphalist, betalist, gammalist, eigVec,
	eigVal))
      printf("matrix2euler: error occured\n");


   if(euler2polar(alphalist,betalist,gammalist,
	thetaGlist,phiGlist,omegaGlist))
      printf("euler2polar: error occured\n");


   if(pca_tensor(eigVal,&da_in,&dr_in))
      printf("pca_tensor: error occured\n");


/*
   if(print_analyzation3d( matrix, eigVec, eigVal, da_in, dr_in, 
	psiList, thetaList, phiList, 
	alphalist, betalist, gammalist,
	thetaGlist, phiGlist, omegaGlist))
      printf("print_analyzation3d: error occured\n");
*/


   *da=da_in;
   *dr=dr_in;

   return 0;
}



/* =================================================================== */
int matrix2euler( matrix, psiList, thetaList, phiList, 
	alphalist, betalist, gammalist, eigenvec,
	eigenvalue)

   Real alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
   Real matrix[9], psiList[8], thetaList[8], phiList[8];
   Real eigenvec[3][3], eigenvalue[3];
{
    int  i, j, n=0, nrot;
    Real sMat[4][4];
    Real eigVec[4][4], eigVal[4], eigVal_temp[4];


    /*** initialization ***/
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

    /* row after row */
    for (j=0; j<3; j++)
       for (i=0; i<3; i++) sMat[i+1][j+1] = matrix[j+3*i];



   /*** determine eigensystem ***/
   jacobimz(sMat,eigVal,eigVec,&nrot);


   /*** sorting (absolute value) of eigenvalues and eigenvectors ***/
   if(sort_eigensystem(eigVal,eigVec))
	printf("sort_eigensystem: error occured\n");


   /*** determine Euler angles according to y-convention ***/
   /*** take into account ambiguity of eigenvectors *(+/-1) ***/
   if(eigVec2euler(1,eigVec,alphalist,betalist,gammalist))
	printf("eigVec2euler: error occured\n");


 

   /*** return values */
   for(i=1;i<=3;i++){
	for(j=1;j<=3;j++)	eigenvec[i-1][j-1] = eigVec[i][j];
	eigenvalue[i-1] = eigVal[i];
   }


    return 0;
}
/* =================================================================== */



/* =============================================================== */
   /*** determine Euler angles according to y-convention ***/
   /*** take into account ambiguity of eigenvectors *(+/-1) ***/
int eigVec2euler(mirror,eigVec,alphalist,betalist,gammalist)

	int mirror;
	Real alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
	Real eigVec[4][4];
{
   int  i,j,n=0;
   Real eigVec_temp[3][3],eigVec_temp2[4][4];
   Real eigVecXlist[NO_LIST][3], eigVecYlist[NO_LIST][3], eigVecZlist[NO_LIST][3];
   Real rot_angle[4];


   /*** shift of matrix ***/
   for(i=0;i<3;i++)
	for(j=0;j<3;j++)	eigVec_temp[i][j]=eigVec[i+1][j+1];


   /*** make list with possible orientations of tensor ***/
   n = make_eigVecList(eigVec_temp,eigVecXlist,eigVecYlist,eigVecZlist);


  /*** analyze the four independent orientations ***/
  for(j=0;j<n;j++) {

  	for(i=0;i<3;i++) {
		eigVec_temp2[i+1][1] = eigVecXlist[j][i];
		eigVec_temp2[i+1][2] = eigVecYlist[j][i];
		eigVec_temp2[i+1][3] = eigVecZlist[j][i];
	}

	if(get_euler(eigVec_temp2,rot_angle))
	   printf("get_euler: error occured\n");

	alphalist[j]=rot_angle[1];
   	betalist[j]=rot_angle[2];
   	gammalist[j]=rot_angle[3];
   }

   return 0;
}
/* =============================================================== */


/* =============================================================== */
int cross_productM(n,m,eigVec,product)

	int n,m;
	Real product[4],eigVec[4][4];
{
   product[1]=eigVec[2][n]*eigVec[3][m]-eigVec[3][n]*eigVec[2][m];
   product[2]=eigVec[3][n]*eigVec[1][m]-eigVec[1][n]*eigVec[3][m];
   product[3]=eigVec[1][n]*eigVec[2][m]-eigVec[2][n]*eigVec[1][m];

   return( 1 );
}
/* =============================================================== */


/* =============================================================== */
int euler2polar(alphalist,betalist,gammalist,
	thetaGlist,phiGlist,omegaGlist)

   Real alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
   Real thetaGlist[NO_LIST], phiGlist[NO_LIST], omegaGlist[NO_LIST];
{
   int i;
   Real a,b,g,temp;


   /* determine the polar angles for the general rotation */
   for(i=0;i<NO_LIST;i++){

	a=RAD(alphalist[i]);
	b=RAD(betalist[i]);
	g=RAD(gammalist[i]);

	phiGlist[i] = PI/2 + (a-g)/2;
	thetaGlist[i] = atan( tan(b/2) / sin((a+g)/2) );
	temp = acos( cos(b/2)*cos((a+g)/2) );

	if( ( sin(b/2)*sin((g-a)/2)/cos(phiGlist[i]) > 0 ) &&
       		( sin(temp)*sin(gammalist[i]) > 0) )
	   	temp = temp;
  	else		
	   	temp = (2*PI - temp);

   	omegaGlist[i] = 2*temp;
   	if(omegaGlist[i]>=2*PI) omegaGlist[i] -= 2*PI;
 
	thetaGlist[i] = DEG(thetaGlist[i]);
	phiGlist[i] = DEG(phiGlist[i]);
	omegaGlist[i] = DEG(omegaGlist[i]);
  }


   return 0;
}
/* =============================================================== */


/* ======== determine magnitude and rhombicity of tensor ========= */
int pca_tensor(eigVal,da,dr)
	Real eigVal[3];
	Real *da, *dr;

{
    Real d[3], rtemp;

    d[0] = eigVal[0]; d[1] = eigVal[1]; d[2] = eigVal[2];
    if (fabs(d[0]) < fabs(d[1])) {RSWAP( d[0], d[1] );}
    if (fabs(d[0]) < fabs(d[2])) {RSWAP( d[0], d[2] );}
    if (fabs(d[1]) < fabs(d[2])) {RSWAP( d[1], d[2] );}

    *da = (d[0] - 0.5*(d[2] + d[1]))/3.0;
    *dr = (d[2] - d[1])/3.0;


   return 0;
}
/* =============================================================== */


/* ===================== print alignment ========================= */
int print_analyzation3d( matrix, eigVec, eigVal, da, dr, 
	psiList, thetaList, phiList, 
	alphalist, betalist, gammalist,
	thetaGlist, phiGlist, omegaGlist)

	Real da, dr;
	Real matrix[9], eigVal[3], eigVec[3][3];
	Real psiList[8],thetaList[8],phiList[8];
	Real alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
	Real thetaGlist[NO_LIST], phiGlist[NO_LIST], omegaGlist[NO_LIST];
{
   int i,j;
   Real sMat[3][3];
   float daHz=0.0,R=0.0;


    /* row after row */
     for (j=0; j<3; j++)
       for (i=0; i<3; i++) sMat[i][j] = matrix[j+3*i];
 

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

/* =============================================================== */
int sort_eigensystem(eigVal,eigVec)
	Real eigVal[4],eigVec[4][4];

{
    int i,j;
    Real rtemp, d[4],vec[4],dvec[4][4];


    /* sort eigenvalues and calculate da and dr */
    for(i=1;i<=3;i++) d[i] = eigVal[i];
    for(i=1;i<=3;i++)
	    for(j=1;j<=3;j++) dvec[i][j] = eigVec[i][j];


    if (fabs(d[1]) < fabs(d[2])) {
	RSWAP( d[1], d[2] );
	for(i=1;i<=3;i++) vec[i]=dvec[i][1];
	for(i=1;i<=3;i++) dvec[i][1]=dvec[i][2];
	for(i=1;i<=3;i++) dvec[i][2]=vec[i];
    }
    if (fabs(d[1]) < fabs(d[3])) {
	RSWAP( d[1], d[3] );
	for(i=1;i<=3;i++) vec[i]=dvec[i][1];
	for(i=1;i<=3;i++) dvec[i][1]=dvec[i][3];
	for(i=1;i<=3;i++) dvec[i][3]=vec[i];
     }
    if (fabs(d[2]) < fabs(d[3])) {
	RSWAP( d[2], d[3] );
	for(i=1;i<=3;i++) vec[i]=dvec[i][2];
	for(i=1;i<=3;i++) dvec[i][2]=dvec[i][3];
	for(i=1;i<=3;i++) dvec[i][3]=vec[i];
    }

    for(i=1;i<=3;i++) eigVal[i] = d[i];
    for(i=1;i<=3;i++)
	    for(j=1;j<=3;j++) eigVec[i][j] = dvec[i][j];


    return 0;
}
/* =============================================================== */

