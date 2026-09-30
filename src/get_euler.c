/*								*/
/* Berechnung der Eulerwinkel aus der Drehmatrix m, d.h. aus	*/
/* Matrix m der normierten Eigenvektoren m=[x,y,z] mit 		*/ 
/* x=x/norm(x), y=y/norm(y) und z=z/norm(z) 			*/
/*								*/
#include <stdlib.h>
#include <math.h>


#define Real double

int get_euler(m,w)
  Real w[4], m[4][4];
{
  int i,j;
  Real norm[4], m_norm[4][4];


   /* normalization */
  for(i=1;i<=3;i++)
  	norm[i] = sqrt(m[1][i]*m[1][i]+m[2][i]*m[2][i]+m[3][i]*m[3][i]);
  for(i=1;i<=3;i++)
	  for(j=1;j<=3;j++)	m_norm[j][i] = m[j][i] / norm[i];


   /* determine euler angles */
   inv_euler(m_norm,w);

   return( 0 );
}
