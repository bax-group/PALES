#include <stdio.h>
#include <string.h>
#include <math.h>


#define NO_LIST 4

#define FPR (void)fprintf
#define RSWAP( RA, RB )  rtemp = RA; RA = RB; RB = rtemp



/* =================================================================== */
int analyzeEigen( eigVal, eigVec, flag )
	int    flag;
	double eigVal[3], eigVec[3][3];
{
    int     i,j,n,error;
    double  alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
    double  eigVecXlist[NO_LIST][3], eigVecYlist[NO_LIST][3];
    double  eigVecZlist[NO_LIST][3];


/***/
/* Initialize:
/***/

    i         = 0;
    j         = 0;
    n         = 0;
    error     = 0;


/***/
/* Sorting of eigenvalues and eigenvectors
/* flag allows selection of sorting criterium
/***/

   if( error = sortEigen(eigVal,eigVec,flag) )
     {
       FPR( stderr, "analyzeEigen: error occured in sortEigen.\n");
       goto shutdown;
     }



/***/
/* Determine Euler angles according to y-convention
/* Take into account ambiguity of eigenvectors *(+/-1)
/***/

   if( eigVec2eulerN(1,eigVec,alphalist,betalist,gammalist,
		     eigVecXlist,eigVecYlist,eigVecZlist) )
     {
        FPR( stderr, "analyzeEigen: error occured in eigVec2euler.\n");
	error = 1;
        goto shutdown;
     }


/***/
/* Output
/***/
    FPR( stderr, "\nanalyzeEigen\n" );
    FPR( stderr, "eigensystem\n" );
    FPR( stderr, "===========\n" );

    FPR( stderr, "eigenvalues (Azz,Ayy,Axx)\n" );
    for(i=0;i<3;i++) FPR( stderr, "%e  ", eigVal[i] );
    FPR( stderr, "\n");
    FPR( stderr, "eigenvectors (z,y,x)\n");
    for(i=0;i<3;i++){
        for(j=0;j<3;j++) FPR( stderr, "%e  ", eigVec[i][j] );
        FPR( stderr, "\n");
    }
    FPR( stderr, "\n" );

    FPR( stderr, "\nOrientation of alignment tensor\n" );
    FPR( stderr, "(4 possible coordinate systems (x/y/z))\n" );
    FPR( stderr, "=======================================\n" );
    FPR( stderr, "x axis  :  ");
    for (i=0;i<NO_LIST;i++) FPR( stderr, "(%.4f/%.4f/%.4f)  ",
				eigVecXlist[i][0],eigVecXlist[i][1],
				eigVecXlist[i][2] );
    FPR( stderr, "\ny axis  :  ");
    for (i=0;i<NO_LIST;i++) FPR( stderr, "(%.4f/%.4f/%.4f)  ",
				eigVecYlist[i][0],eigVecYlist[i][1],
				eigVecYlist[i][2] );
    FPR( stderr, "\nz axis  :  ");
    for (i=0;i<NO_LIST;i++) FPR( stderr, "(%.4f/%.4f/%.4f)  ",
				eigVecZlist[i][0],eigVecZlist[i][1],
				eigVecZlist[i][2] );
    FPR( stderr, "\n\n");
    FPR( stderr, "(Euler angles [deg.] for 4 equivalent orientations)\n" );
    FPR( stderr, "==================================================\n" );
    FPR( stderr, "alpha   :  ");
    for (i=0;i<NO_LIST;i++) FPR( stderr, "%7.2f  ",alphalist[i] );
    FPR( stderr, "\nbeta    :  ");
    for (i=0;i<NO_LIST;i++) FPR( stderr, "%7.2f  ",betalist[i] );
    FPR( stderr, "\ngamma   :  ");
    for (i=0;i<NO_LIST;i++) FPR( stderr, "%7.2f  ",gammalist[i] );
    FPR( stderr, "\n\n");


/***/
/* Exit point:
/***/

shutdown:

    return( error );
}
/* =================================================================== */


/* =============================================================== */
int sortEigen(eigVal,eigVec,flag)
	int    flag;
	double eigVal[3],eigVec[3][3];

{
    int     i,j,error;
    double  d[3],vec[3],dvec[3][3],rtemp;


/***/
/* Initialize:
/***/

    i         = 0;
    j         = 0;
    error     = 0;


/***/
/* Sort eigenvalues according to flag
/***/
    for(i=0;i<3;i++) d[i] = eigVal[i];
    for(i=0;i<3;i++)
	    for(j=0;j<3;j++) dvec[i][j] = eigVec[i][j];


  if( flag )
  {

    if (fabs(d[0]) < fabs(d[1])) {
	RSWAP( d[0], d[1] );
	for(i=0;i<3;i++) vec[i]=dvec[i][0];
	for(i=0;i<3;i++) dvec[i][0]=dvec[i][1];
	for(i=0;i<3;i++) dvec[i][1]=vec[i];
    }
    if (fabs(d[0]) < fabs(d[2])) {
	RSWAP( d[0], d[2] );
	for(i=0;i<3;i++) vec[i]=dvec[i][0];
	for(i=0;i<3;i++) dvec[i][0]=dvec[i][2];
	for(i=0;i<3;i++) dvec[i][2]=vec[i];
     }
    if (fabs(d[1]) < fabs(d[2])) {
	RSWAP( d[1], d[2] );
	for(i=0;i<3;i++) vec[i]=dvec[i][1];
	for(i=0;i<3;i++) dvec[i][1]=dvec[i][2];
	for(i=0;i<3;i++) dvec[i][2]=vec[i];
    }

  }

  else {


    if (d[0] < d[1]) {
	RSWAP( d[0], d[1] );
	for(i=0;i<3;i++) vec[i]=dvec[i][0];
	for(i=0;i<3;i++) dvec[i][0]=dvec[i][1];
	for(i=0;i<3;i++) dvec[i][1]=vec[i];
    }
    if (d[0] < d[2]) {
	RSWAP( d[0], d[2] );
	for(i=0;i<3;i++) vec[i]=dvec[i][0];
	for(i=0;i<3;i++) dvec[i][0]=dvec[i][2];
	for(i=0;i<3;i++) dvec[i][2]=vec[i];
     }
    if (d[1] < d[2]) {
	RSWAP( d[1], d[2] );
	for(i=0;i<3;i++) vec[i]=dvec[i][1];
	for(i=0;i<3;i++) dvec[i][1]=dvec[i][2];
	for(i=0;i<3;i++) dvec[i][2]=vec[i];
    }

  }

    for(i=0;i<3;i++) eigVal[i] = d[i];
    for(i=0;i<3;i++)
	    for(j=0;j<3;j++) eigVec[i][j] = dvec[i][j];


/***/
/* Exit point:
/***/

shutdown:

    return( error );
}
/* =============================================================== */


/* =============================================================== */
   /*** determine Euler angles according to y-convention ***/
   /*** take into account ambiguity of eigenvectors *(+/-1) ***/
int eigVec2eulerN(mirror,eigVec,alphalist,betalist,gammalist,
		  eigVecXlist,eigVecYlist,eigVecZlist)

   int    mirror;
   double alphalist[NO_LIST], betalist[NO_LIST], gammalist[NO_LIST];
   double eigVec[3][3];
   double eigVecXlist[NO_LIST][3], eigVecYlist[NO_LIST][3], eigVecZlist[NO_LIST][3];
{
   int    i,j,n=0,error;
   double eigVec_temp[3][3],eigVec_temp2[4][4];
   double rot_angle[4];



/***/
/* Initialize:
/***/

    i         = 0;
    j         = 0;
    n         = 0;
    error     = 0;

    for(i=0;i<3;i++)
	for(j=0;j<3;j++)	eigVec_temp[i][j]=eigVec[i][j];



/***/
/* Make list with possible orientations of tensor
/***/

   n = make_eigVecList(eigVec_temp,eigVecXlist,eigVecYlist,eigVecZlist);
   if( n != NO_LIST )
     {
       FPR( stderr, "eigVec2eulerN: error occured in make_eigVecList.\n" );
       error = 1;
       goto shutdown;
     }


/***/
/* Get Euler angles for different coordinate systems
/***/

  for(j=0;j<n;j++) {

  	for(i=0;i<3;i++) {
		eigVec_temp2[i+1][1] = eigVecXlist[j][i];
		eigVec_temp2[i+1][2] = eigVecYlist[j][i];
		eigVec_temp2[i+1][3] = eigVecZlist[j][i];
	}

	if( get_euler(eigVec_temp2,rot_angle) )
	{
           FPR( stderr, "eigVec2eulerN: error occured in get_euler.\n" );
	   error = 1;
           goto shutdown;
	}

	alphalist[j] = rot_angle[1];
   	betalist[j]  = rot_angle[2];
   	gammalist[j] = rot_angle[3];
   }


/***/
/* Exit point:
/***/

shutdown:

    return( error );
}
/* =============================================================== */


