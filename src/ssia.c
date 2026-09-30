#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "pales.h"
#include "memory.h"
#include "prec.h"

#include "hsymbol.h"

#define RAD(x) (x*PI/180.)
#define DEG(x) (x*180./PI)
#define FPR (void)fprintf

#define epsilon		7.82853e+01  /* unitless */
#define k_B		8.31451e-03  /* [kJ/mol.K] */
#define e_2		1.38935e+02  /* [nm.kJ/mol] */

double Szz(double tensor_zaxis[4]);
double Sxx(double tensor_xaxis[4]);
double Syy(double tensor_yaxis[4]);
double Sxy(double tensor_xaxis[4], double tensor_yaxis[4]);
double Sxz(double tensor_xaxis[4], double tensor_zaxis[4]);
double Syz(double tensor_yaxis[4], double tensor_zaxis[4]);



int stericAlignment( palesInfo, rimOrder, lcType, 
        stepOutLength, max_obstruction,
        radius,gridSpacing,wvAmount,molarRatio,
	no_atoms,x,y,z,
        max_protein_length,min_protein_length,
        digit_psi, n_dots, dots, saupe )

	int    lcType, no_atoms, digit_psi, n_dots;
        double rimOrder, stepOutLength;
        double radius,gridSpacing,wvAmount,molarRatio;
	double max_obstruction, *x, *y, *z, *dots;
        double max_protein_length, min_protein_length;
        double saupe[5];
        struct PalesInfo *palesInfo;
{
  int     i,j,k,add_pts;
  int     total_bic=0,total_ph=0;
  int     no_ok_pts_bic=0,no_ok_pts_ph=0;
  double  atom_vec[4],atom_vec_temp[4];
  double  radius2=0.0;
  double  scale_fac=0.0;
  static  double saupe_dis[MAX_NO_PT][5];
  double  saupe_bic[5],saupe_ph[5];
  int     m=0, orientation_flag=1;

  double  x_axis_tensor[4],y_axis_tensor[4],z_axis_tensor[4];
  double  Dreh_m[4][4];
  float   *pot_dist;
  double  *total;

  int     a,t_pts,no_pt_ph,lengthGrid;
  int     no_add_orient, no_add_ph_orient;
  float   offsetGrid,*grid_el;
  double  d,p,end_dist;


  int     kk,mm;
  double  inc_psi;
  double  beta,gamma,alpha;

  double  szz, sxx_yy, sxy, sxz, syz;
  double  tensor_xaxis[4], tensor_yaxis[4], tensor_zaxis[4];
  int     o_pts;

  double  saupe_dev[5];
  double  sum1,sum2,sum3,sum4,sum5,sum6;


  int     error;
  int     bic_phage;
  double  density;
  float   *xJ,*yJ;


/***/
/* Initialize:
/***/

    error       = 0;
    density     = 1.46;
    bic_phage   = 0;

    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    sum6 = 0.0;
    offsetGrid = 0;

    xJ         = (float *)NULL;
    yJ         = (float *)NULL;
    pot_dist   = (float *)NULL;
    grid_el    = (float *)NULL;
    total      = (double *)NULL;


/***/
/* Allocate and extract coord lists, etc:
/***/

    xJ       = fltAlloc( "stericAlignment",  no_atoms+1 );
    yJ       = fltAlloc( "stericAlignment",  no_atoms+1 );
    pot_dist = fltAlloc( "stericAlignment",  MAX_NO_PT+1 );
    grid_el  = fltAlloc( "stericAlignment",  MAX_NO_PT );
    total    = dFltAlloc( "stericAlignment", MAX_NO_PT );

    if (!xJ || !yJ || !pot_dist || !grid_el || !total)
       {
           FPR( stderr, "Error allocating coord buffers.\n" );
           error = 1;
           goto shutdown;
       }


/***/
/* Determine deviation from uniform sampling of orientations:
/***/

   inc_psi = 360.0/digit_psi;
   for( i=0; i<n_dots; i++ )
      {

        /* mapping from three-dimensional Cartesian coordinates */
        /* to spherical coordinates as in cart2sph.m */
        alpha = atan2(dots[i*3+1],dots[i*3]);
        beta  = atan2(dots[i*3+2],sqrt(dots[i*3]*dots[i*3]+dots[i*3+1]*dots[i*3+1]));
        beta  = PI/2.0-beta;

        for( mm=0; mm<digit_psi; mm++ )
           {                   
              gamma = RAD(0.0+mm*inc_psi);

              rotate_coorsystem( beta,alpha,gamma,tensor_xaxis,
                        tensor_yaxis,tensor_zaxis );

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

  o_pts = n_dots*digit_psi;
  saupe_dev[0] = sum1/(1.0*o_pts);
  saupe_dev[1] = (sum2-sum3)/(1.0*o_pts);
  saupe_dev[2] = sum4/(1.0*o_pts);
  saupe_dev[3] = sum5/(1.0*o_pts);
  saupe_dev[4] = sum6/(1.0*o_pts);

/*
    FPR( stderr, "Saupe dev matrix\n" );
    for(i=0;i<5;i++) FPR( stderr, "%e  ",saupe_dev[i] );
    FPR( stderr, "\n" );
*/


  /*************/
  /*** setup initial tranlational grid ***/
   for( i=0; i<MAX_NO_PT; i++ ) grid_el[i]=0.0;

   offsetGrid = min_protein_length + max_obstruction;

   t_pts = (max_protein_length-min_protein_length)/gridSpacing;
   if( t_pts >= MAX_NO_PT )
     {
       t_pts =  MAX_NO_PT - 1;
       FPR( stderr, "WARNING: concentration too low - magnitude of alignment tensor incorrect\n" );
     }

   for( i=0; i<=t_pts; i++ )
      {
        grid_el[i] = offsetGrid + i * gridSpacing;
      }
   no_pt_ph = t_pts;



  /*************/
  /*** calculate scaling ***/
   if( lcType == 1 )
     {
        no_add_orient = (stepOutLength-grid_el[t_pts])/gridSpacing;
     }
   if( lcType==2 )
     {
        end_dist = max_obstruction/sqrt(wvAmount/density); 
        a= (end_dist-max_obstruction) / gridSpacing;
        no_add_ph_orient = a + (gridSpacing/max_obstruction) * (a*(a+1)/2.0);
     }



  /*************/
  /*** plane ***/
   if( lcType == 1 )
     {

       /*** initialization ***/
        for( j=0; j<=t_pts; j++ )
           {
              pot_dist[j] = radius-(grid_el[j]-max_obstruction);
           }

        for( i=0; i<5; i++ )
                for(j=0;j<=t_pts;j++)
                        saupe_dis[j][i] = 0.0;
        for( j=0; j<=t_pts; j++ )
                total[j] = 0;



        /*** perform simulation ***/
        kk = 0;
        for(i=0; i<n_dots; i++){

           /* mapping from three-dimensional Cartesian coordinates */
           /* to spherical coordinates as in cart2sph.m */
           alpha = atan2(dots[i*3+1],dots[i*3]);
           beta = atan2(dots[i*3+2],sqrt(dots[i*3]*dots[i*3]+dots[i*3+1]*dots[i*3+1]));
           beta = PI/2.0-beta;

           for(mm=0;mm<digit_psi;mm++) {                        
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

                        xJ[k]=0;
                        for(j=1;j<=3;j++)       xJ[k]+=Dreh_m[j][1]*atom_vec[j];
                }
                

                /* sort x coordinates into ascending order */
                /* slows the program down */
                /* shell(no_atoms,xJ); */


                /* check hard potential for different distances */
                for(j=0;j<=t_pts;j++) {

                        orientation_flag=1;
                        for(m=1;m<=no_atoms;m++) {
                                if(xJ[m] <= pot_dist[j]) { 
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
   if( ((lcType==1) && (bic_phage==1)) || (lcType==2) )
     {
        
      /* initialization */
        radius2 = max_obstruction + radius;


        /* set up complete grid */
        end_dist = max_obstruction/sqrt(wvAmount/density)-max_obstruction; 
        i = 0;
        d = grid_el[no_pt_ph];
        while( d<=end_dist )
             {  
                i++;
                grid_el[no_pt_ph+i] = grid_el[no_pt_ph] + gridSpacing*i;
                d= grid_el[no_pt_ph+i];

                if( no_pt_ph+i == (MAX_NO_PT-1) )
                  {
                      FPR( stderr, "WARNING: concentration too low - magnitude of alignment tensor incorrect\n" );
                      break;
                  }

             }
        lengthGrid = i+no_pt_ph;
        /* lengthGrid = no_pt_ph; */


        for(i=0;i<5;i++)
                for(j=0;j<=lengthGrid;j++)
                        saupe_dis[j][i] = 0.0;
        for(j=0;j<=lengthGrid;j++)
                total[j]=0.0;


        /*** perform simulation ***/
        kk = 0;
        for( i=0; i<n_dots; i++ )
           {

           /* mapping from three-dimensional Cartesian coordinates */
           /* to spherical coordinates as in cart2sph.m */
           alpha = atan2(dots[i*3+1],dots[i*3]);
           beta = atan2(dots[i*3+2],sqrt(dots[i*3]*dots[i*3]+dots[i*3+1]*dots[i*3+1]));
           beta = PI/2.0-beta;

           for(mm=0;mm<digit_psi;mm++) {                        
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

                        xJ[k]=0;
                        yJ[k]=0;
                        for(j=1;j<=3;j++){
                                xJ[k]+=Dreh_m[j][1]*atom_vec[j];
                                yJ[k]+=Dreh_m[j][2]*atom_vec[j];
                        }
                }


                /* shift center of mass of protein to new position on grid */
                for(j=0;j<=lengthGrid;j++) {

                        p=1.0;
                        orientation_flag=1;
                        if(j<=no_pt_ph)
                        for(m=1;m<=no_atoms;m++) {
                                if( sqrt((xJ[m]+grid_el[j])*(xJ[m]+grid_el[j])+(yJ[m]*yJ[m])) <= radius2) { 
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
           {
                for(j=0;j<=lengthGrid;j++) {
                        if(total[j]>0.0)
                                saupe_dis[j][i] = saupe_dis[j][i] / total[j] - saupe_dev[i];
                        else
                                saupe_dis[j][i] = 0.0;
                }
           }



        no_ok_pts_ph = 0;
        total_ph     = 0;
        for(j=0;j<=lengthGrid;j++)
           {
                if( total[j] != 0.0 )       no_ok_pts_ph++;
                total_ph += total[j];
           }

        /* for(i=0;i<=lengthGrid;i++) printf("grid_el[%d]: %.1f  total[%d]: %.1f     saupe_dis[%d][0]: %e\n",
                                        i,grid_el[i],i,total[i],i,(1.0/(o_pts*1.0)) *saupe_dis[i][0]); */

        /* calculate Saupe matrix for range of grid points */
        /* for which simulation is performed */
        for( i=0; i<5; i++ )
           {
                saupe_ph[i] = 0;
                for( j=0; j<=lengthGrid; j++)
                   {
                        d = (total[j]/(o_pts*1.0)) * saupe_dis[j][i] * grid_el[j]/max_obstruction;
                        /* if(i==0) printf("grid_el[%d]/max_obstruction: %.1f   saupe_dis[%d][0]: %e   saupe_disSCALE[%d][0]: %e\n",
                                j,grid_el[j]/max_obstruction,j,saupe_dis[j][0],j,d); */
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
    if ((lcType==1) && (bic_phage==1.0))
        for(i=0;i<5;i++) {
                saupe[i] = 1.0/(1+(1.0/molarRatio)) * (saupe_bic[i] -  
                (1.0/molarRatio)*((1.0*no_ok_pts_bic) / (1.0*no_ok_pts_ph)) * rimOrder * saupe_ph[i]);
        }
    if ((lcType==1) && (bic_phage!=1.0))
        for(i=0;i<5;i++) {
                saupe[i] = saupe_bic[i];
                /* saupe[i] = 1.0/(1+(1.0/molarRatio)) * (1-(1.0/molarRatio) * rimOrderrimOrder) * saupe_bic[i];
*/
        }
    if( lcType == 2 )
      {
        for( i=0; i<5; i++ )
           {
                saupe[i] = saupe_ph[i];
           }
      }


    /* normalize and scale saupe matrix */
    /* allow for possibility not to scale down the saupe matrix */
    scale_fac = 1.0;

    if ( lcType==2 )
       {
         add_pts = no_add_ph_orient;
       }
    else
       {
         add_pts = total_bic/o_pts + no_add_orient;
       }


    for( i=0; i<5; i++ )
       {
         saupe[i] = saupe[i] / (add_pts*scale_fac);
       }

/*
    FPR( stderr, "Saupe matrix\n" );
    for(i=0;i<5;i++) FPR( stderr, "%e  ",saupe[i] );
    FPR( stderr, "\n" );
*/


/***/
/* Exit point:
/***/

shutdown:

    (void) deAlloc( "stericAlignment", xJ,       sizeof(float)*(no_atoms+1) );
    (void) deAlloc( "stericAlignment", yJ,       sizeof(float)*(no_atoms+1) );
    (void) deAlloc( "stericAlignment", pot_dist, sizeof(float)*(MAX_NO_PT+1) );
    (void) deAlloc( "stericAlignment", grid_el,  sizeof(float)*(MAX_NO_PT) );
    (void) deAlloc( "stericAlignment", total,    sizeof(double)*(MAX_NO_PT) );


    return( error );
}
/* ========================================================================= */


/* ================ simulate macromolecular alignment ========================= */
int macroAlignment( palesInfo, rimOrder, lcType, 
        stepOutLength, max_obstruction,
        radius,gridSpacing,wvAmount,molarRatio, temperatureKPot,
        chargeInfo,kPot,kEner,saltConc,surfChDens,cutOffPOT,
	no_atoms,x,y,z,
        digit_psi, n_dots, ddots,
        max_protein_length,min_protein_length,
        chCount, chargeAtomID, charges,
        digPotential, potential, offset, stepsize,
        saupe )

	int    lcType, no_atoms, digit_psi, n_dots;
        int    digPotential, chCount, *chargeAtomID;
        float  *charges;
        double rimOrder, stepOutLength;
        double radius,gridSpacing,wvAmount,molarRatio;
	double max_obstruction, *x, *y, *z, *ddots;
        double max_protein_length, min_protein_length;
        double offset,stepsize,*potential;
        double saupe[5];
        struct PalesInfo *palesInfo;
{
  int     i,j,k,error;
  double  atom_vec[4],atom_vec_temp[4];
  double  radius2=0.0;
  static  double  saupe_dis[MAX_NO_PT][5];
  int     m=0, orientation_flag=1;

  double  x_axis_tensor[4],y_axis_tensor[4],z_axis_tensor[4];
  double  Dreh_m[4][4];
  double  cosd[4];
  double  *total;
  double  total_bic=0.0, total_ph=0.0;
  double  endPotential, p, potvalue, pot_energy;
  float   *xJ,*yJ,*grid_el;
  int     t_pts, no_pt_ph;
  double  d, dist;
  int     lengthGrid, id, rID;
  double  end_dist,density;

  int     kk,mm;
  float   inc_psi;
  double  beta,gamma,alpha,test;

  int     o_pts;
  double  szz, sxx_yy, sxy, sxz, syz;
  double  tensor_xaxis[4], tensor_yaxis[4], tensor_zaxis[4];

  double  saupe_dev[5];
  double  sum1,sum2,sum3,sum4,sum5,sum6;
  float   cutoff, offsetGrid;

  double  rho,kappa,surf_pot,tan_surf_pot,avo;


/***/
/* Initialize:
/***/

    error       = 0;
    density     = 1.46;
    avo         = 6.022*1e+23;
    sum1        = 0.0;
    sum2        = 0.0;
    sum3        = 0.0;
    sum4        = 0.0;
    sum5        = 0.0;
    sum6        = 0.0;
    offsetGrid  = 0;
    cutoff      = 100.0;

    xJ         = (float *)NULL;
    yJ         = (float *)NULL;
    grid_el    = (float *)NULL;
    total      = (double *)NULL;


/***/
/* Allocate and extract coord lists, etc:
/***/

    xJ       = fltAlloc( "macroAlignment",  no_atoms+1 );
    yJ       = fltAlloc( "macroAlignment",  no_atoms+1 );
    grid_el  = fltAlloc( "macroAlignment",  MAX_NO_PT );
    total    = dFltAlloc( "macroAlignment", MAX_NO_PT );

    if (!xJ || !yJ || !grid_el || !total)
       {
           FPR( stderr, "Error allocating coord buffers.\n" );
           error = 1;
           goto shutdown;
       }


/***/
/* Determine deviation from uniform sampling of orientations:
/***/
FPR( stderr, "@@2: ddots[0] = %e\n", ddots[0] );

  o_pts   = n_dots * digit_psi;
  inc_psi = 360.0  / digit_psi;
/*
  for( i=0; i<n_dots; i++ )
     {

        alpha = atan2( ddots[i*3+1], ddots[i*3] );
        beta  = atan2( ddots[i*3+2], sqrt(ddots[i*3]*ddots[i*3]+ddots[i*3+1]*ddots[i*3+1]));
        beta  = PI/2.0 - beta;

        for( mm=0; mm<digit_psi; mm++ )
           {                   
                gamma = RAD(0.0+mm*inc_psi);

                rotate_coorsystem(beta,alpha,gamma,tensor_xaxis,
                        tensor_yaxis,tensor_zaxis);

                sum1 += Szz(tensor_zaxis);
                sum2 += Sxx(tensor_xaxis);
                sum3 += Syy(tensor_yaxis);
                sum4 += Sxy(tensor_xaxis,tensor_yaxis);
                sum5 += Sxz(tensor_xaxis,tensor_zaxis);
                sum6 += Syz(tensor_yaxis,tensor_zaxis);
          }
    }
  saupe_dev[0] = sum1/(1.0*o_pts);
  saupe_dev[1] = (sum2-sum3)/(1.0*o_pts);
  saupe_dev[2] = sum4/(1.0*o_pts);
  saupe_dev[3] = sum5/(1.0*o_pts);
  saupe_dev[4] = sum6/(1.0*o_pts);
  for( i=0; i<5; i++ )
        FPR( stderr, "o_pts: %d    saupe_dev[%d]: %e\n",
               o_pts,i,saupe_dev[i] );
*/


  /*************/
  /*** setup initial tranlational grid ***/
   for( i=0; i<MAX_NO_PT; i++ ) grid_el[i]=0.0;

   offsetGrid = min_protein_length + max_obstruction;
   t_pts = (max_protein_length-min_protein_length)/gridSpacing;
   for( i=0; i<=t_pts; i++)
        grid_el[i] = offsetGrid + i * gridSpacing;
   no_pt_ph=t_pts;



  /*************/
  /*** plane ***/
  if( lcType == 1 ) {

        /*** initialization ***/
        for( j=0; j<=t_pts; j++)   radius2 = radius + max_obstruction;


        /* set up complete grid */
        if( stepOutLength - cutOffPOT > 0 )
                beta = cutOffPOT;
        else
                beta = stepOutLength;

        i=0;
        d=grid_el[t_pts];
        while( d <= beta )
             { 
               i++;
               if( t_pts+i == MAX_NO_PT )
                  {
                      FPR( stderr, "WARNING: concentration too low - magnitude of alignment tensor incorrect\n" );
                      break;
                  }
 
               grid_el[t_pts+i] = grid_el[t_pts] + gridSpacing*i;
               d = grid_el[t_pts+i];
            }
        lengthGrid = i+t_pts;
        /* printf("lengthGrid: %d\n",lengthGrid); */


        for( i=0; i<5; i++ )
                for(j=0;j<=lengthGrid;j++)
                        saupe_dis[j][i] = 0.0;
        for( i=0; i<=lengthGrid; i++ )
           {
              total[i]=0.0;
           }


        /*** calculate electrostatic parameters ***/
        rho      = avo * 1.0e-24 * saltConc;
        kappa    = sqrt(8*PI*rho*e_2/(epsilon*k_B*temperatureKPot));
        test     = surfChDens * kappa / (4.0*rho);
        surf_pot = 2.0 * log(test+sqrt(test*test+1.0));
        tan_surf_pot = tanh(surf_pot/4.0);
        if( chargeInfo )
          {
            printf("\ncharacteristics of electrostatic potential:\n");
            printf("===========================================\n");
            printf("equilibrium density of ion pairs [1/nm_3]: %e\n",rho);
            printf("Debye length [nm]: %e\n",1./kappa);
            printf("potential at surface [e/kB*T]: %e\n\n\n",surf_pot);
          }


        /*** perform simulation ***/
        kk = 0;
        for( i=0; i<n_dots; i++ )
           {

           /* mapping from three-dimensional Cartesian coordinates */
           /* to spherical coordinates as in cart2sph.m */
           alpha = atan2( ddots[i*3+1], ddots[i*3] );
           beta  = atan2( ddots[i*3+2], sqrt(ddots[i*3]*ddots[i*3]+ddots[i*3+1]*ddots[i*3+1]));
           beta  = PI/2.0 - beta;

           for(mm=0;mm<digit_psi;mm++) {                        
                gamma = RAD(0.0+mm*inc_psi);

                /* report how many orientations have been tested */
                /* if(fmod(++kk,1000)==0.0)     printf("# of orientations: %d\n",kk); */
                /* if(fmod(++kk,200)==0.0)      printf(".\n");  */


                /* calculate elements of Saupe matrix for Euler angles i, */
                /* e.g. polar angles theta (beta) and phi (alpha) */
                rotate_coorsystem(beta,alpha,gamma,tensor_xaxis,
                        tensor_yaxis,tensor_zaxis);

                szz    = Szz(tensor_zaxis);
                sxx_yy = Sxx(tensor_xaxis)-Syy(tensor_yaxis);
                sxy    = Sxy(tensor_xaxis,tensor_yaxis);
                sxz    = Sxz(tensor_xaxis,tensor_zaxis);
                syz    = Syz(tensor_yaxis,tensor_zaxis);


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

                        xJ[k]=0;
                        for(j=1;j<=3;j++)
                                xJ[k]+=Dreh_m[j][1]*atom_vec[j];
                }


                /* shift center of mass of protein to new position on grid */
                for( j=0; j<=lengthGrid; j++) {

                        /*** steric evaluation ***/
                        orientation_flag=1;
                        if(j<=t_pts) {
                           for(m=1;m<=no_atoms;m++) {
                                if((xJ[m]+grid_el[j]) <= radius2) { 
                                        orientation_flag=0;
                                        p=0.0;
                                        break;
                                }
                           }
                        }


                        /*** electrostatic evaluation ***/
                        if( orientation_flag ) {
                                pot_energy = 0.0;
                                for(m=1;m<=chCount;m++) {
                                        id = chargeAtomID[m];
                                        dist = (xJ[id]+grid_el[j])-max_obstruction;
                                        /* test = -0.9*exp(-1.*dist/20) */
                                        test = tan_surf_pot*exp(-1.*kappa*dist/10.0);
                                        potvalue = 4.0* 0.5*log( (1.+test)/(1.-test) );
/*
                                        rID = dist/stepsize+1;
                                        potvalue = ((potential[rID+1]-potential[rID])/stepsize)*(dist-stepsize*(rID-1))+potential[rID];
*/
                                        pot_energy += charges[m] * potvalue;
                                }
                                p = exp((-1.0)*kEner*pot_energy);
                        }

                        /*** calculate order matrix (not normalized) ***/                       
                        saupe_dis[j][0] += p*szz;
                        saupe_dis[j][1] += p*sxx_yy;
                        saupe_dis[j][2] += p*sxy;
                        saupe_dis[j][3] += p*sxz;
                        saupe_dis[j][4] += p*syz;

                        total[j]+=p;    

                } /* all translational positions */

           } /* end gamma loop */
        } /* end dots loop */


        total_bic = 0.0;
        for( j=0; j<=lengthGrid; j++) total_bic += total[j];
        if( stepOutLength-cutOffPOT > 0 )
                total_bic += o_pts*(stepOutLength-cutOffPOT)/gridSpacing;

        /* for(i=0;i<=lengthGrid;i++) printf("grid_el[%d]: %.1f  total[%d]: %.1f     saupe_dis[%d][0]: %e\n",
                                        i,grid_el[i],i,total[i],i,saupe_dis[i][0]); */
        


        /* sum up over translational grid points */
        for( i=0; i<5; i++)
           {
                saupe[i] = 0.0;
                for(j=0;j<=lengthGrid;j++)
                        saupe[i] += saupe_dis[j][i];
                saupe[i] /= total_bic;
           }

   }  /*** end plane part ***/




  /****************/
  /*** cylinder ***/
  if( lcType == 2 ) {


        /*** initialization ***/
        radius2 = max_obstruction + radius;

        /*************/
        /*** refine electrostatic potential ***/
        for( i=1; i<=digPotential; i++ ) potential[i] *= kPot;




        /* set up complete grid */
        end_dist = max_obstruction/sqrt(wvAmount/density)-max_obstruction; 
        endPotential = (digPotential-1)*stepsize+offset;
        if (end_dist<(endPotential+max_protein_length))
                printf("warning: potential of neighbouring liquid crystal particles overlapping!");

        if(end_dist-cutOffPOT>0)
                beta = cutOffPOT;
        else
                beta = end_dist;


        i=0;
        d=grid_el[no_pt_ph];
        while( d <= beta )
             {  
                i++;

                if( no_pt_ph+i == MAX_NO_PT )
                  {
                      FPR( stderr, "WARNING: concentration too low - magnitude of alignment tensor incorrect\n" );
                      break;
                  }

                grid_el[no_pt_ph+i] = grid_el[no_pt_ph] + gridSpacing*i;
                d= grid_el[no_pt_ph+i];
             }
        lengthGrid = i+no_pt_ph;
        /* printf("lengthGrid: %d\n",lengthGrid); */


        for(i=0;i<5;i++)
                for(j=0;j<=lengthGrid;j++)
                        saupe_dis[j][i] = 0.0;
        for(i=0;i<=lengthGrid;i++)      total[i]=0.0;



        /*** perform simulation ***/
        kk=0;
        for(i=0; i<n_dots; i++){

           /* mapping from three-dimensional Cartesian coordinates */
           /* to spherical coordinates as in cart2sph.m */
           alpha = atan2(ddots[i*3+1],ddots[i*3]);
           beta = atan2(ddots[i*3+2],sqrt(ddots[i*3]*ddots[i*3]+ddots[i*3+1]*ddots[i*3+1]));
           beta = PI/2.0-beta;

           for(mm=0;mm<digit_psi;mm++) {                        
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

                        xJ[k]=0;
                        yJ[k]=0;
                        for(j=1;j<=3;j++){
                                xJ[k]+=Dreh_m[j][1]*atom_vec[j];
                                yJ[k]+=Dreh_m[j][2]*atom_vec[j];
                        } 

                        /* rotate protein instead of coordinate system */
                        /* rotate_vector(beta,alpha,gamma,atom_vec,atom_vec_temp);
                        xJ[k]=atom_vec_temp[1];
                        yJ[k]=atom_vec_temp[2]; */
                }

                /* shift center of mass of protein to new position on grid */
                for(j=0;j<=lengthGrid;j++) {

                        /*** steric evaluation ***/
                        orientation_flag=1;
                        if(j<=no_pt_ph) {
                           for(m=1;m<=no_atoms;m++) {
                                if( sqrt((xJ[m]+grid_el[j])*(xJ[m]+grid_el[j])+(yJ[m]*yJ[m])) <= radius2 ) { 
                                        orientation_flag=0;
                                        p=0.0;
                                        break;
                                }
                           }
                        }


                        /*** electrostatic evaluation ***/
                        if(orientation_flag) {
                                pot_energy=0.0;
                                for(m=1;m<=chCount;m++) {
                                        id = chargeAtomID[m];
                                        dist = sqrt((xJ[id]+grid_el[j])*(xJ[id]+grid_el[j])+yJ[id]*yJ[id])-max_obstruction;

                                        /* wall potential */
                                        /* test = -0.9*exp(-dist/20.0);
                                        potvalue = 1.8* 0.5*log( (1.+test)/(1.-test) ); */

                                        /* cylinder potential */
                                        rID = dist/stepsize+1;
                                        /* linear interpolation */
                                        potvalue = ((potential[rID+1]-potential[rID])/stepsize)*(dist-stepsize*(rID-1))+potential[rID];
                                        pot_energy += charges[m] * potvalue;
                                }
                                p = exp((-1.0)*kEner*pot_energy);
                        }


                        /*** calculate order matrix (not normalized) ***/                       
                        saupe_dis[j][0] += p*szz;
                        saupe_dis[j][1] += p*sxx_yy;
                        saupe_dis[j][2] += p*sxy;
                        saupe_dis[j][3] += p*sxz;
                        saupe_dis[j][4] += p*syz;

                        total[j]+=p;    

                } /* all translational positions */

           } /* end gamma loop */
        } /* end dots loop */




        total_ph=0.0;
        for(j=0;j<=lengthGrid;j++) total_ph+=total[j]*grid_el[j]/max_obstruction;
        if(end_dist>cutOffPOT) {
                i=0;
                d=grid_el[lengthGrid];
                while(d<=end_dist) {  
                        i++;
                        d = grid_el[lengthGrid] + gridSpacing*i;
                        total_ph+=o_pts*d/max_obstruction;
                }
        }

        /* for(i=0;i<=lengthGrid;i++) printf("grid_el[%d]: %.1f  total[%d]: %.1f     saupe_dis[%d][0]: %e\n",
                                        i,grid_el[i],i,total[i],i,(1.0/total_ph) *saupe_dis[i][0]);
        */


        /* sum up over translational grid points */
        /* '*grid_el[j]/max_obstruction' taking into account */
        /* the increasing number of allowed states with */
        /* increasing distance from the phage surface */
        for(i=0;i<5;i++) {
                saupe[i] = 0.0;
                for(j=0;j<=lengthGrid;j++)
                        saupe[i] += (1.0/total_ph)*saupe_dis[j][i] * grid_el[j]/max_obstruction;
        }

    } /* end phage part */



   /*** output ***/
    /* printf("Saupe\n");
    for(i=0;i<5;i++) printf("%e  ",saupe[i]);
    printf("\n"); */
 


/***/
/* Exit point:
/***/

shutdown:

    (void) deAlloc( "macroAlignment", xJ,       sizeof(float)*(no_atoms+1) );
    (void) deAlloc( "macroAlignment", yJ,       sizeof(float)*(no_atoms+1) );
    (void) deAlloc( "macroAlignment", grid_el,  sizeof(float)*(MAX_NO_PT) );
    (void) deAlloc( "macroAlignment", total,    sizeof(double)*(MAX_NO_PT) );


    return( error );
}
/* ========================================================================= */


/* ======= calculation of Szz according to Saupe ===================== */
double Szz(double tensor_zaxis[4])
{
        return 0.5 * (3 * tensor_zaxis[3]*tensor_zaxis[3] - 1.0);
}
/* ========================================================================= */


/* ========== calculation of Sxx according to Saupe =============== */
double Sxx(double tensor_xaxis[4])
{
        return 0.5 * (3 * tensor_xaxis[3]*tensor_xaxis[3] - 1.0);
}
/* ========================================================================= */


/* ========== calculation of Syy according to Saupe =============== */
double Syy(double tensor_yaxis[4])
{
        return 0.5 * (3 * tensor_yaxis[3]*tensor_yaxis[3] - 1.0);
}
/* ========================================================================= */


/* ============= calculation of Sxy according to Saupe ============== */
double Sxy(double tensor_xaxis[4], double tensor_yaxis[4])
{
        return 0.5 * 3 * tensor_xaxis[3] * tensor_yaxis[3];
}
/* ========================================================================= */


/* ============ calculation of Sxz according to Saupe =============== */
double Sxz(double tensor_xaxis[4], double tensor_zaxis[4])
{
        return 0.5 * 3 * tensor_xaxis[3] * tensor_zaxis[3];
}
/* ========================================================================= */


/* ============ calculation of Syz according to Saupe =============== */
double Syz(double tensor_yaxis[4], double tensor_zaxis[4])
{
        return 0.5 * 3 * tensor_yaxis[3] * tensor_zaxis[3];
}
/* ========================================================================= */
