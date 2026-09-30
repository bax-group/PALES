/***/
/* dc.h: Dipolar coupling and Chemical shift homology and simulation.
/***/

#define DC_FIRSTTIME_YES          1
#define DC_FIRSTTIME_NO           0

#define DC_LIST_YES               1
#define DC_LIST_NO                0
#define DC_UVAL               -9999
#define DC_UNDEF               9999

#define DC_UPDATE_YES             1
#define DC_UPDATE_NO              0

#define DC_SVD                    1
#define DC_SAUPE                  2
#define DC_FIT_NONLIN             3
#define DC_GRIDSEARCH             4
#define DC_FIXED_ALL              5
#define DC_STERIC                 6
#define DC_ELECTROSTATIC          7

/* PALES 2.1 names for the two variable-tensor fitting modes. */
#define DC_FIXED_NONE             3
#define DC_FIXED_DADR             4

/* PALES selection modes (not present in newer NMRPipe headers). */
#define DC_EXCLUDE                1
#define DC_INCLUDE_ALL            2
#define DC_INCLUDE_MAX            3

#define MFR_DC_SVD               1
#define MFR_DC_FIXED_DA          2
#define MFR_DC_FIXED_RH          3
#define MFR_DC_FIXED_DADR        4
#define MFR_DC_FIXED_DARH        5
#define MFR_DC_FIXED_DADR_ROT    6

#define DC_MAXTENSORS            32 
#define PCS_MAXTENSORS           32 

#define CS_UNDEF_TEST( VAL )    ((VAL) >  900)
#define DC_UNDEF_TEST( VAL )    ((VAL) >  900)

#define DC_MIN_DIST 0.01

/***/
/* Methyl tetrahedral angle theta = 110.9
/* DC methyl motional averaging scale is (3*cos2(theta) - 1)/2.
/***/

#define DC_METHYL_AVG -0.309107

/***/
/* Physical constants:
/***/

#define hCons     6.625581e-4
#define hbCons    1.05457266e-4

#define gH        2.675198e+8
#define gN       -2.7116e+7
#define gC        6.7283e+7
#define gD        4.1066e+7
#define gF        2.5170e+8
#define gP        1.0829e+8
#define gFE       8.680624e+06 

#define RAD_N_HN      1.041
#define RAD_CO_HN     2.055
#define RAD_CO_N      1.329
#define RAD_CO_CA     1.525
#define RAD_CA_HA     1.107
#define RAD_N_CA      1.458
#define RAD_CH        1.107
#define RAD_CH3       1.117
#define RAD_CH_ARO    1.097

#define DC_CLASS_UNDEF  -1
#define DC_CLASS_NONE    0
#define DC_CLASS_HN_N    1
#define DC_CLASS_HN_HA   2
#define DC_CLASS_HN_C    3
#define DC_CLASS_HN_CA   4
#define DC_CLASS_N_CA    5
#define DC_CLASS_N_C     6
#define DC_CLASS_HA_CA   7
#define DC_CLASS_HA_HB   8
#define DC_CLASS_CA_CB   9
#define DC_CLASS_C_CA   10 
#define DC_CLASS_CH3    11 

#define HMS_NONE         0
#define HMS_ADDED        1
#define HMS_SUBTRACTED  -1

#define DIVAL_HN      -21585.20

/***/
/* Dipolar Coupling Info:
/***/

struct DInfo 
   {
    int    id;
    int    absFlag;
    int    count;
    int    valid;
    int    selected;
    int    class;
    int    resIDI;
    int    resIDJ;
    char   resNameI[PDB_STRLEN];
    char   resNameJ[PDB_STRLEN];
    char   atomNameI[PDB_STRLEN];
    char   atomNameJ[PDB_STRLEN];
    char   chainNameI[PDB_STRLEN];
    char   chainNameJ[PDB_STRLEN];
    char   segNameI[PDB_STRLEN];
    char   segNameJ[PDB_STRLEN];
    float  d;
    float  dCalc;
    float  ddc;
    float  ddcCalc;
    float  di;
    float  diCalc;
    float  w;
    float  wCalc;
    float  sigma;
    float  sigmaCalc;
    float  s;
    float  fc;
    float  s1, s2, s3, s4, s5;
    int    nVec;
    int    diFixedFlag;
    struct AtomInfo *atomI[DYN_MAXDEGENERACY]; 
    struct AtomInfo *atomJ[DYN_MAXDEGENERACY]; 
    int    atomICount;
    int    atomJCount;
   };

struct DTensor
   {
    float daRef, drRef, rhRef;
    float psiRef, thetaRef, phiRef; 
    float psiRel, thetaRel, phiRel; 
    float psiOff, thetaOff, phiOff; 
    float psiAdd, thetaAdd, phiAdd; 
    float da, daMin, daMax;
    float dr, drMin, drMax;
    float rh, rhMin, rhMax;
    float gridStep;
    float psiList[8],   psi,   psiMin,   psiMax;
    float thetaList[8], theta, thetaMin, thetaMax;
    float phiList[8],   phi,   phiMin,   phiMax;
    float saupe[5], dSaupe[5], userSaupe[5];
    float eigVal[3];
    float sXYZ[3];
    float dcRMS;
    float fc;
    int   eulerFlag;
    int   eulerCount;
    int   tensorMode;
    int   rotRefFlag;
    int   autoFlag;
    int   autoSignFlag;
    int   rhFlag;
    int   daScaleFlag;
    int   bestIter;
    Matrix4x4 m;
   };

struct DTable
   {
    struct DTensor dTensor;

    /* PALES 2.1 keeps observed and calculated values in parallel arrays. */
    struct DInfo   *dObsList;
    struct DInfo   *dCalcList;
    int    dCount;

    /* Newer NMRPipe spelling, retained for source compatibility. */
    struct DInfo   *dcList;
    struct PDBInfo **pdbList;
    int    dcCount;
    int    pdbCount;
    int    validCount;
    int    svdCount;
    int    dMinRes;
    int    dMaxRes;
    int    nullChainName;
    int    nullSegName;
    int    maxCount;
    int    workSizeR;
    int    workSizeI;
    float  *workR;
    int    *workI;
   };

struct DTableList
   {
    struct DTable *dTab[DC_MAXTENSORS];
    int    dTabCount;
   };

/***/
/* Chemical Shift Info:
/***/

struct CSInfo
   {
    int    id;
    int    count;
    int    valid;
    int    selected;
    int    predictable;
    int    resID;
    char   resName[PDB_STRLEN];
    char   resName1[PDB_STRLEN];
    char   resName3[PDB_STRLEN];
    char   atomName[PDB_STRLEN];
    char   chainName[PDB_STRLEN];
    char   segName[PDB_STRLEN];
    float  w;
    float  x;
    float  y;
    float  z;
    float  fc;
    float  shift;
    float  predShift;
    float  secShift;
    float  surfShift;
    float  rcShift;
    float  hmShift;
    float  sigma;
    float  gamma;
    struct AtomInfo      *atomList[DYN_MAXDEGENERACY];
    int    atomCount;
    struct ResInfo       *res;
    struct PhiPsiSurface *surf;
    char   *auxInfo;
   };

struct PCSInfo
   {
    int    id;
    int    valid;
    int    selected;
    int    resID;
    char   resName[PDB_STRLEN];
    char   atomName[PDB_STRLEN];
    char   chainName[PDB_STRLEN];
    char   segName[PDB_STRLEN];
    float  x;
    float  y;
    float  z;
    float  fc;
    struct AtomInfo *atom;
    float  gamma;
    struct DTensor dTensor;
   };

struct CSTable
   {
    struct SeqInfo *seqInfo;
    struct CSInfo  *csList;
    struct PCSInfo *pcsList;
    int    csCount;
    int    pcsCount;
    int    validCSCount;
    int    validPCSCount;
    int    csMinRes;
    int    csMaxRes;
    int    nullChainName;
    int    nullSegName;
    int    nullPCSChainName;
    int    nullPCSSegName;
    int    workSize;
    int    hmsMode;
    float  *work;
   };

struct PhiPsiSurface *findCSSurface();

/* Stand-alone PALES compatibility interface. */
int dTensorNull(struct DTensor *tensor);
int dTabNull(struct DTable *table);
int dTabFree(struct DTable *table);
int csTabNull(struct CSTable *table);
int csTabFree(struct CSTable *table);
int readDC(char *fileName, struct DTable *table);
int mapDCInfo(struct PDBInfo *pdbInfo, struct DTable *table);
int selectD(struct DInfo *list, int count, int first, int last, int mode);
int buildDC(struct PDBInfo *pdbInfo, struct DTable *table, int offset,
            int segmentFirst, float *weights, int weightCount);
int fitD(struct PDBInfo *pdbInfo, struct DTable *table, int verbose);
int getStatDC(struct DInfo *observed, struct DInfo *calculated, int count,
              float *rms, float *chi2, int *validCount);
int rot2euler(float rotation[3][3], float *psiList,
              float *thetaList, float *phiList);
int euler2saupe(float *saupe, float psi, float theta, float phi,
                float da, float dr);
int reportEuler(FILE *stream, struct DTensor *tensor);
