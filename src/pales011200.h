#define MAX_M 400			/*  # of rows, max # of dipolar couplings  */
#define MAX_NO_DC 2000
#define	MAX_ATOM_NO 5000
#define	MAX_NO_PT 5000
#define	MAX_NO_PT_SQR 5000
#define	MAX_ORIENT_PT 1000
#define	RADIUS 1.0
#define MAXCHAR  100				/* maximale Laenge der Pfadnamen */
#define DC_THRESH  1000				/* cutoff for allowed DCs */
#define MAX_DIG_POT 10000
#define MAX_DOTS 3000
#define MAX_ENSEMBLE 40000
#define MaxCountCharge 800


/***/
/* PALES modes:
/***/

#define PALES_STERIC              1
#define PALES_STERIC_FREE         2
#define PALES_STERIC_GRID         3
#define PALES_ELECTROSTATIC       4
#define PALES_ELECTROSTATIC_FREE  5
#define PALES_ELECTROSTATIC_GRID  6
#define PALES_DC                  7
#define PALES_CS                  8
#define PALES_HIST                9
#define PALES_ASSIGN              10
#define PALES_DC_CP               11
#define PALES_SAUPE_CP            12


#define	MAXCHAR  100
#define MAX_OPTIONS 10

#define FABS(X) (((X) < 0.0) ? -(X) : (X))

#define NO_LIST 4
#define H_NO 329
#define N_ind 211
#define HN_ind 551
#define CA_ind 114
#define HA_ind 476
#define C_ind 1


int debug=0;			/* debugging switched off */


struct ElInfo
   {
    int    elFlag;              /* flag for incorporation of electrostatics (0/1)*/
    int    pkaFlag;		/* using pKa default values */
    float  pH;			/* original pH_auto */
    float  concNaCl;		/* [M] molarity mol/l */
    float  chDensSurf;		/* surface charge density [e/nm2] */
    int    chTerminiFlag;	/* flag for putting charges at termini (0/1)*/
    double kEnergy;		/* scaling for EL energy */
    double kPot;		/* scaling for EL potential */
    float  kNegCh;		/* scaling down of negative charges */
    float  kPosCh;		/* scaling down of positive charges */
    float  cutoffPot;		/* cutoff for taking into account */
				/* any interaction (Angstrom) */
    float  temperaturePot;	/* temperature [K] */
    int    chFlag;
   };

struct GridInfo
   {
    int	   outEflag;		/* put out surface ensemble */
    int	   densE;        	/* density of dots on ensemble sphere */
    int    angleDomFlag;		/* change rel. orientation of domains (0/1) */
    float  angleDom[2];         /* Euler angles for adjusting */
					     /* angle between 2 domains */
    int    resDom1Origin; 		/* origin for orientation of domain 1 */
    int	   resDom2Origin; 		/* origin for orientation of domain 2 */
    int	   startResDom2,startResDom1,endResDom1,endResDom2;
					/* domain boundaries */
    int    digPsiE; 		/* digitization of rotation about */
					/* 3rd Euler angle */
    float  lBetaCone;
    float  hBetaCone;
    float  cAlphaCone;
    float  errAlphaCone;
    int    startResPart1,endResPart1,startResPart2,endResPart2;
    int    clashFlag;		/* check steric clash of domains (0/1) */
    float  clashDist;		/* distance for steric clash */
    int    initOrientFlag;
    int    pivotCenter;		/* 1: translate molecule to center of mass */
    int	   fixRotDomain;		/* 1: rotate protein such that domain which */
   };

struct PalesInfo
   {
    int    lcType;
    double concentration;
    int    densityNSC;
    int    digPsiNSC;
    double radiusModel;
    double gridSpacing;
    double includeRim, molarRatio, lengthModel, orderDiscRim;
    double lcOrder;		/* LIQUID_CRYSTALL_ORDER: [0-1] */
    int    selFlag;		/* complete protein is taken into account (1/2/12/3) */
    int    Hflag;			/* flag for taking protons into account (0/1) */
    double radiusAtom;		/* radius of pseudo atoms [A] */
    double rSingleCell;
    struct GridInfo gridInfo;
    struct ElInfo   elInfo;
   };


/* @@ */
float   Pmax=10.0;
/* @@ */

int    liquid_type=1;		/* bicelle */
double bic_phage=0;		/* take rim of bicelle into account */
int    grid_type=2;		/* 1D grid */
double concentration=0.05;	/* total lipid concentration [%] */
double ellipticity=1.0;		/* not used */
double max_obs_dim=20.0;	/* radius of bicelle rim or of phage (30.0) [A] */
double grid_spacing=0.5;	/* spacing of grid [A] */
double medium_dist_max=8000.0;	/* not used */
double molar_ratio=3.0;		/* DMPC/DHPC molar ratio */
double order_parameter=0.8;	/* LIQUID_CRYSTALL_ORDER: [0-1] */
int    select_flag=0;		/* complete protein is taken into account (1/2/12/3) */
int    refine_ok=1,scale_down_ok=1; /* REFINE_SAUPE: & SCALE_SAUPE: */
int    mode=0;			/* newest program mode (0/1/2)*/
int    H_flag=0;		/* flag for taking protons into account (0/1) */
double atom_radius=3.0;		/* radius of pseudo atoms [A] */
double radius_bic=200, order_rel=0.25;
double step_out_length;


/*** output ***/
int svd_flag=1;			/* flag for doing SVD */
int charge_info=1;


/*** electrostatics ***/
int electrostatics=0;		/* flag for incorporation of electrostatics (0/1)*/
int pKa_default=0;		/* using pKa default values */
int ch_termini=1;		/* flag for putting charges at termini (0/1)*/
double k_ener=1.0;		/* scaling for EL energy */
double k_pot=1.0;		/* scaling for EL potential */
float ch_weight=1.0;		/* weight charges */
float ch_neg_scale=1.0;		/* scaling down of negative charges */
float ch_pos_scale=1.0;		/* scaling down of positive charges */
float pH_auto=5.4;
float cutoffPOT = 100.0;	/* cutoff for taking into account */
				/* any interaction (Angstrom) */
float salt_conc = 0.05;		/* [M] molarity mol/l */
float surf_ch_dens = 0.8;	/* surface charge density [e/nm2] */
float temperaturePot =  298.150;	/* temperature [K] */


/*** digitization of spherical grid ***/
int old_mode = 0;
int dot_density = 100;
int dig_theta=9,dig_phi=18,dig_psi=18;


/*** ensemble averaging ***/
int	ensemble_out=0;		/* put out surface ensemble */
int	ensemble_averaging=0;	/* ensemle averaging flag (0/1) */
int	ensemble_density=200;	/* density of dots on ensemble sphere */
int     change_angle=0;		/* change rel. orientation of domains (0/1) */
float   alphaDom2=0,betaDom2=0,gammaDom2=0;  /* Euler angles for adjusting */
					     /* angle between 2 domains */
int     resDom1Origin=0; 	/* origin for orientation of domain 1 */
int	resDom2Origin=0; 	/* origin for orientation of domain 2 */
int	startResDom2,startResDom1,endResDom1,endResDom2;
				/* domain boundaries */
int     digit_ensemble=1; 	/* digitization of rotation about */
				/* 3rd Euler angle */
/* limits for polar angles of cone */
float	betaConeL = 0.0;
float   betaConeH = 30.0;
float   alphaConeC = 0.0;
float   alphaConeErr = 180.0;
float   alphaConeCInv = 0.0;

int     startResPart1,endResPart1,startResPart2,endResPart2;

int     checkMode = 1;		/* check steric clash of domains (0/1) */
float   clashDistance = 3.6;	/* distance for steric clash */
int     initOrientMode = 1;
int     pivotCenter = 1;	/* 1: translate molecule to center of mass */
				/* prior to PALES prediction */
int	fixRotDomain = 0;	/* 1: rotate protein such that domain which */
				/* is rotated is returned to starting position */
