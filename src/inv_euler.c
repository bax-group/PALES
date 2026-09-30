/*								*/
/* Berechnung der Eulerwinkel aus der Drehmatrix m, d.h. aus	*/
/* Matrix m der normierten Eigenvektoren m=[x,y,z] mit 		*/ 
/* x=x/norm(x), y=y/norm(y) und z=z/norm(z) 			*/
/*								*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#define pi 3.14159265

#define Real double



void inv_euler(m,w)

  Real w[4], m[4][4];
{
  int i;
  Real beta, s;
  Real cos_alpha,sin_alpha,alpha_grad;
  Real cos_gamma,sin_gamma;
  Real cos_alphgam,sin_alphgam,gamma_grad;


  for(i=1;i<=3;i++)	w[i]=0.0;

  beta=acos(m[3][3]);

/*--------------------------- beta > 0 -------------------------*/
  if( fabs(sin(beta)) > 1.0e-10 ) {

    cos_alpha = m[3][1]/sin(beta);
    sin_alpha = m[3][2]/sin(beta);
    if(sin_alpha >= 0) {
      if(cos_alpha >= 0) 
         alpha_grad = acos(cos_alpha)*180/pi;           /* 1. Quadrant */
      else
         alpha_grad  = acos(cos_alpha)*180/pi;          /* 2. Quadrant */
    }
    else {
      if(cos_alpha <= 0.0)
        alpha_grad = 360 - acos(cos_alpha)*180/pi;      /* 3. Quadrant */
      else
        alpha_grad = 360 - acos(cos_alpha)*180/pi;      /* 4. Quadrant */
    }

  /* jetzt dito mit gamma */
    cos_gamma  = -m[1][3] / sin(beta);
    sin_gamma  =  m[2][3] / sin(beta);
    if(sin_gamma >= 0){
      if(cos_gamma >= 0)
        gamma_grad = acos(cos_gamma)*180/pi;
      else
        gamma_grad = acos(cos_gamma)*180/pi;
    }
    else {
      if(cos_gamma <= 0.0)   
        gamma_grad = 360 - acos(cos_gamma)*180/pi;   
      else
        gamma_grad = 360 - acos(cos_gamma)*180/pi;
    }

  }
  else {
    /* ----- für beta < 1e-12  */
      s=1; 
      if(m[3][3] < 0)	s=-1;
      cos_alphgam = s*m[1][1];    /* alphgam steht entweder für alpha + gamma oder alpha - gamma */
      sin_alphgam = s*m[1][2];
      if(sin_alphgam >= 0) {
        if(cos_alphgam >= 0)   
          alpha_grad = acos(cos_alphgam)*180/pi;
        else
          alpha_grad = acos(cos_alphgam)*180/pi;
      }
      else {
        if(cos_alphgam <= 0)
          alpha_grad = 360  - acos(cos_alphgam)*180/pi;   
        else
          alpha_grad = 360  - acos(cos_alphgam)*180/pi;
      }

      gamma_grad = 0;
  }

  w[1]=alpha_grad;
  w[2]=beta*180/pi;
  w[3]=gamma_grad;

}
