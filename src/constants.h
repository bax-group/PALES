#define RAD(x) (x*PI/180.)
#define DEG(x) (x*180./PI)
#define S(x) sin(x)
#define C(x) cos(x)
#define PI (3.141592654)

#define gh (2.675198e4)
#define gn ((-1.0)*(2.7116e3))
#define gc (6.7283e3)
#define h (1.054494e-27)
#define csa (-160.0)

/* #define rnhn (1.042)	*/	/* from Ottiger and Bax, JACS 120, 12334 */
#define rnhn (1.042)
#define rcaha (1.119)
#define rnco (1.329)
#define rcohn (2.069)
#define rcaco (1.525)
#define rhh (2.3)		/* randomly chosen */
#define facnhn (gn*gh/(rnhn*rnhn*rnhn))
#define normcaha ( facnhn / (gc*gh/(rcaha*rcaha*rcaha)) )
#define normnco ( facnhn / (gn*gc/(rnco*rnco*rnco)) )
#define normcohn ( facnhn / (gc*gh/(rcohn*rcohn*rcohn)) )
#define normcaco ( facnhn / (gc*gc/(rcaco*rcaco*rcaco)) )
#define normhh ( facnhn / (gh*gh/(rhh*rhh*rhh)) )


/* #define hCons     6.625581e-4
#define hbCons    1.054494e-4 */
#define hCons     6.6260755e-4		/* as given in Bronstein */
#define hbCons    1.05457266e-4		/* as given in Bronstein */
#define u0	  PI*4.0e-7
#define cx	  ((u0*hbCons)/(4.0*PI*PI))

#define gH        2.675198e+8
#define gN       -2.7116e+7
#define gC        6.7283e+7
#define gD        4.1066e+7
#define gP        1.0841e+8

#define RAD_N_HN   1.020	/* does not incorporate fast angular motions */
#define RAD_CO_HN  2.041
#define RAD_CO_N   1.329
#define RAD_CO_CA  1.525
#define RAD_CA_HA  1.095
#define RAD_N_CA   1.458

double dcNORMFACTORS[6] = {1.0,normcaha,normcaco,normnco,normcohn,normhh};


#define epsilon		7.82853e+01  /* unitless */
#define k_B		8.31451e-03  /* [kJ/mol.K] */
#define e_2		1.38935e+02  /* [nm.kJ/mol] */
