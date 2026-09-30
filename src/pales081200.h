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

#define OP_ADD 1
#define OP_SUB 2
#define OP_MUL 3
#define OP_DIV 4
#define OP_MAX 5
#define OP_MIN 6
#define OP_EMB 7
#define OP_IL  8
#define OP_MAC 9


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


struct ElInfo
   {
    int    elFlag;              /* flag for incorporation of electrostatics (0/1)*/
    int    pkaFlag;		/* using pKa default values */
    float  pH;			/* original pH_auto */
    float  concNaCl;		/* [M] molarity mol/l */
    float  chDensSurf;		/* surface charge density [e/nm2] */
    int    chTerminiFlag;	/* flag for putting charges at termini (0/1)*/
    float  kEnergy;		/* scaling for EL energy */
    float  kPot;		/* scaling for EL potential */
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
    int    palesMode;
    int    lcType;
    float  concentration;
    int    densityNSC;
    int    digPsiNSC;
    float  radiusModel;
    float  gridSpacing;
    float  includeRim, molarRatio, lengthModel, orderDiscRim;
    float  lcOrder;		/* LIQUID_CRYSTALL_ORDER: [0-1] */
    int    Hflag;			/* flag for taking protons into account (0/1) */
    int    surfFlag;
    float  radiusAtom;		/* radius of pseudo atoms [A] */
    float  rSingleCell;
    int    opCode;
    float  c1, c2;
    int    eulerFlag;
    int    verbFlag;
    struct GridInfo gridInfo;
    struct ElInfo   elInfo;
   };

struct DcInfo
   {
    int    segFirst;
    int    segLast;
    int    segID;
    int    segStep;
    int    aFirst; 
    int    excludeMode;
   };
