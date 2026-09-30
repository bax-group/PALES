/***/
/* pdb.h: some definitions for manipulating PDB data.
/***/

#define PDB_STRLEN                16
#define PDB_NAME_MAX              32

#define INIT_ATOM_COUNT           4096
#define INCR_ATOM_COUNT           1024

#define INIT_RES_COUNT            512 
#define INCR_RES_COUNT            150

#define INIT_RESATOM_COUNT        32 
#define INCR_RESATOM_COUNT        8 

#define INIT_CHAIN_COUNT          8
#define INCR_CHAIN_COUNT          8 

#define INIT_SEG_COUNT            8
#define INCR_SEG_COUNT            8 

#define PDB_HTAB_SIZE             4111 

#define PDB_BAD_RESID_VAL      -666666

#define PDB_UNDEF_LIM         666666.0
#define PDB_UNDEF_VAL         777777.0

#define PDB_SELECT_NULL       (struct PDBSelectInfo *)NULL

#define PDB_MAP_NO             0
#define PDB_MAP_YES            1

#define PDB_NAMING_OLD         0
#define PDB_NAMING_IUPAC       1

#define PDB_ALIGN_NONE         0
#define PDB_ALIGN_ALL          1
#define PDB_ALIGN_TRANSLATE    2

#define PDB_SELECT_INV        -1
#define PDB_SELECT_CLEAR       0
#define PDB_SELECT_NONE        0
#define PDB_SELECT_ALL         1
#define PDB_SELECT_SET         1
#define PDB_SELECT_PARSE       2
#define PDB_SELECT_PROTEIN     4
#define PDB_SELECT_REGION      5
#define PDB_SELECT_SPHERE      6

#define PDB_EXCLUDE            1
#define PDB_INCLUDE_ALL        2
#define PDB_INCLUDE_MIN        3
#define PDB_INCLUDE_MAX        4

#define SURFACE_NONE           0
#define SURFACE_SHIFT          1
#define SURFACE_CONFORMATION   2

#define SURFACE_PSIPREV       -1
#define SURFACE_PHIPSI         0
#define SURFACE_PHINEXT        1

#define DYN_MAXCHI            30
#define DYN_MAXRADGYR         32
#define DYN_MAXDEGENERACY      9
#define DYN_MAXSHIFTTYPES      7

#define PDB_BAD_RESID_TEST( ID ) ((ID) <= PDB_BAD_RESID_VAL)
#define PDB_UNDEF_TEST( VAL )    ((VAL) > PDB_UNDEF_LIM)

struct Coord3D
   {
    float x;
    float y;
    float z;
   };

struct DCoord3D
   {
    double x;
    double y;
    double z;
   };

struct PhiPsiSurface
   {
    char  surfName[PDB_STRLEN];
    float *value;
    float *density;
    float *dev;
    float minVal;
    float maxVal;
    float avgVal;
    float avgDev;
    float thresh;
    int   type;
    int   offsetType;
    int   obsCount;
    int   phiCount;
    int   psiCount;
    float phi1;
    float phiN;
    float psi1;
    float psiN;
   };

#define DYN_NPREV 3

struct AtomInfo
   {
    int    atomID;
    char   *auxInfo;
    char   atomName[PDB_STRLEN];
    char   locID[PDB_STRLEN];
    char   resName[PDB_STRLEN];
    char   chainID[PDB_STRLEN];
    int    resID;
    char   icode[PDB_STRLEN];
    float  x, y, z;
    float  occupancy;
    float  tempFactor;
    char   segID[PDB_STRLEN];
    char   element[PDB_STRLEN];
    char   chargeStr[PDB_STRLEN];
    int    model;
    int    selected;
    struct ResInfo   *resInfo;
    struct ChainInfo *chainInfo;
    struct SegInfo   *segInfo;
    struct DCoord3D  velocity[DYN_NPREV];
    struct DCoord3D  force[DYN_NPREV];
    double time[DYN_NPREV];
    float  ljEpsilon;
    float  ljSigma;
    float  vdwfc;
    float  radius;
    float  charge;
    float  mass;
    float  shift;
    float  dshift;
    float  rcShift;
    int    fixed;
    int    hetFlag;
   };

struct ResInfo
   {
    int   resID;
    int   atomCount;
    int   maxAtomCount;
    char  *auxInfo;
    char  resName[PDB_STRLEN];
    char  chainID[PDB_STRLEN];
    char  segID[PDB_STRLEN];
    struct ChainInfo *chainInfo;
    struct SegInfo   *segInfo;
    struct ResInfo   *next;
    struct ResInfo   *prev;
    float phi;
    float psi;
    float val;
    struct ResTypeInfo   *resType;
    struct PhiPsiSurface *phiPsiSurf;
    struct AtomInfo      **atomList;
   };

#ifdef MAC_OSX
struct ResTypeInfo
   {
    char *resName;
    char resLetter;
    int  resCode;
    int  resClass;
   };
#else
struct ResTypeInfo
   {
    char resName[PDB_STRLEN+1];
    char resLetter;
    int  resCode;
    int  resClass;
   };
#endif

struct ChainInfo
   {
    int   atomCount;
    int   resCount;
    int   maxAtomCount;
    int   maxResCount;
    int   minResID;
    int   maxResID;
    char  *auxInfo;
    float val;
    char  chainID[PDB_STRLEN];
    char  segID[PDB_STRLEN];
    struct AtomInfo **atomList;
    struct ResInfo  **resList;
    struct ResInfo  **resArray;
    struct SegInfo   *segInfo;
   };

struct SegInfo
   {
    int   atomCount;
    int   resCount;
    int   minResID;
    int   maxResID;
    int   maxAtomCount;
    int   maxResCount;
    char  *auxInfo;
    float val;
    char  segID[PDB_STRLEN];
    char  chainID[PDB_STRLEN];
    struct AtomInfo   **atomList;
    struct ResInfo    **resList;
    struct ResInfo    **resArray;
    struct ChainInfo  *chainInfo;
   };

struct PDBInfo
   {
    int    pdbID;
    char   *auxInfo;
    struct AtomInfo  *atomList;
    struct ResInfo   *resList;
    struct ChainInfo *chainList;
    struct SegInfo   *segList;
    struct Coord3D   *coordList;
    struct Coord3D   *vdwList;
    HTable *atomTab;
    HTable *resTab;
    HTable *chainTab;
    HTable *segTab;
    int    atomCount;
    int    resCount;
    int    chainCount;
    int    segCount;
    int    modelCount;
    int    minResID;
    int    maxResID;
    int    maxAtomCount;
    int    maxResCount;
    int    maxChainCount;
    int    maxSegCount;
    int    maxCoordCount;
    int    firstAtomLine; 
    int    lastAtomLine; 
    char   **bondMatrix;
    int    *bondCount;
    int    **neighborList;
    int    *neighborCount;
    char   defChainName[PDB_STRLEN+1];
    char   defSegName[PDB_STRLEN+1];
   };

struct PDBSelectInfo
   {
    int   res1, resN;
    char  resNameList[NAMELEN+1]; 
    char  atomNameList[NAMELEN+1]; 
    char  chainNameList[NAMELEN+1]; 
    char  segNameList[NAMELEN+1]; 
    float x1, xn; 
    float y1, yn; 
    float z1, zn; 
    float rad;
    int   resNameCount;
    int   atomNameCount;
    int   chainNameCount;
    int   segNameCount;
    int   resIDMode;
    int   resNameMode;
    int   atomNameMode;
    int   chainNameMode;
    int   segNameMode;
    int   selectID;
    int   selectMode;
    int   selectClass;
    int   entry1;
    int   entryN;
    int   useRecordTypeFlag;
    int   hetatmFlag;
    int   atomFlag;
    int   waterFlag;
    int   dbgFlag;
   };

struct pdbFieldInfo
   {
    int  fieldLoc;
    char *fieldName;
    int  loc1;
    int  locn; 
   };

struct PDBSeqInfo 
   {
    char pdbName[NAMELEN+1];
    char *seq;
    int  r1, rN;
    int  seqOff;
   };

struct GenericEnergy
   {
    int    type;
    int    status;
    double energy;
    double delta;
    double diff;
    double val;
    double vio;
    double aux;
   };

/***/
/* PDB Fields:
/***/

#define ATOMID_LOC           0
#define ATOMNAME_LOC         1
#define LOCID_LOC            2
#define RESNAME_LOC          3
#define CHAINID_LOC          4
#define RESID_LOC            5
#define ICODE_LOC            6
#define X_LOC                7
#define Y_LOC                8
#define Z_LOC                9
#define OCCUPANCY_LOC       10
#define TEMPFACTOR_LOC      11
#define SEGID_LOC           12
#define ELEMENT_LOC         13
#define CHARGE_LOC          14

#ifdef USE_PDB_FIELDS
static struct pdbFieldInfo pdbFieldList[] =
   {
    ATOMID_LOC,     "ATOMID",       7, 11, 
    ATOMNAME_LOC,   "ATOMNAME",    13, 16,
    LOCID_LOC,      "LOCID",       17, 17,
    RESNAME_LOC,    "RESNAME",     18, 21,
    CHAINID_LOC,    "CHAINID",     22, 22,
    RESID_LOC,      "RESID",       23, 26,
    ICODE_LOC,      "ICODE",       27, 27,
    X_LOC,          "X",           31, 38,
    Y_LOC,          "Y",           39, 46,
    Z_LOC,          "Z",           47, 54,
    OCCUPANCY_LOC,  "OCCUPANCY",   55, 60,
    TEMPFACTOR_LOC, "TEMPFACTOR",  61, 66,
    SEGID_LOC,      "SEGID",       73, 76,
    ELEMENT_LOC,    "ELEMENT",     77, 78,
    CHARGE_LOC,     "CHARGE",      79, 80,
    0,              (char *)NULL,  0,  0, 
   };
#endif

int  pdbReset(), indexPDB(), pdbInit();
int  resNull(), chainNull(), segNull();

/* Stand-alone PALES interface. */
int readPDB(char *fileName, struct PDBInfo *pdbInfo,
            struct PDBSelectInfo *selection);
int writePDB(char *fileName, struct PDBInfo *pdbInfo,
             struct PDBSelectInfo *selection, char *remark);
int pdbFree(struct PDBInfo *pdbInfo);
int pdbNull(struct PDBInfo *pdbInfo);
int atomNull(struct AtomInfo *atomInfo);
int pdbSelectNull(struct PDBSelectInfo *selection);
int pdbSetSelection(struct PDBInfo *pdbInfo, int selected);
int pdbSelect(struct PDBInfo *pdbInfo, struct PDBSelectInfo *selection);
int pdbSelect1(struct PDBInfo *pdbInfo, int res1, int resN,
               char *atomNames, int selectID, int matchMode);
int getPDBSelectArgs(int argc, char **argv,
                     struct PDBSelectInfo *selection);
int storePDBCoords(struct PDBInfo *pdbInfo);
int restorePDBCoords(struct PDBInfo *pdbInfo);

struct AtomInfo  **getNextAtomPtr();
struct AtomInfo  *getNextAtom();
struct ResInfo   *getNextRes();
struct ChainInfo *getNextChain();
struct SegInfo   *getNextSeg();

struct AtomInfo  *findAtom(), *findAtomN();
struct ResInfo   *findRes();
struct ChainInfo *findChain();
struct SegInfo   *findSeg();

struct PDBSeqInfo *pdbSeqAlloc();

struct ResTypeInfo *findResType(), *findResTypeByName();
char   *findResName(), *getResName1(), *getResName3();

double dynDist(), dynAngle(), dynDihedral();
float  pdbDist(), pdbAngle(), pdbDihedral();

float  getPDBVDWRadius(), getPDBAtomMass();

int    updatePDBAngles(struct PDBInfo *pdbInfo);
int    updatePDBShifts(struct PDBInfo *pdbInfo);
char   *pdbStrCpy(), *pdbNameCpy();

#include "pdb3d.h"

