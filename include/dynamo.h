/******************************************************************************/
/*                                                                            */
/*              ---- DYNAMO Molecular Dynamic Software System ----            */
/*                       Copyright 1999 Frank Delaglio                        */
/*                                                                            */
/*               This software is not for distribution without                */
/*                  the written permission of the author.                     */
/*                                                                            */
/******************************************************************************/

/***/
/* dynamo.h: definitions for DYNAMO.
/***/

#define DYNAMO_REV_TXT "2013.047.09.17"
#define DYNAMO_VER_TXT "7.2"

#define DYN_CHECK_KEY       1
#define DYN_IGNORE_KEY      2

#define DYN_READ_RAMA       1
#define DYN_IGNORE_RAMA     2

#define DYN_STATUS_VIO     -1
#define DYN_STATUS_NONE     0
#define DYN_STATUS_OK       1
#define DYN_STATUS_UNUSED   2
#define DYN_STATUS_INVALID  3

#define DYN_TYPE_AUTO      -1 
#define DYN_TYPE_NONE       0
#define DYN_TYPE_PDB        1
#define DYN_TYPE_GMC        2
#define DYN_TYPE_GDB        3
#define DYN_TYPE_ATOM       4
#define DYN_TYPE_DC         5
#define DYN_TYPE_PCS        6
#define DYN_TYPE_CS         7
#define DYN_TYPE_J          8
#define DYN_TYPE_DIST       9
#define DYN_TYPE_NOE       10 
#define DYN_TYPE_HMS       11 
#define DYN_TYPE_HBOND     12 
#define DYN_TYPE_DS        13 
#define DYN_TYPE_AC        14 
#define DYN_TYPE_BOND      15
#define DYN_TYPE_ANGLE     16  
#define DYN_TYPE_UANGLE    17  
#define DYN_TYPE_TORSION   18
#define DYN_TYPE_IMPROPER  19 
#define DYN_TYPE_VDW_EXCL  20 
#define DYN_TYPE_VDW       21 
#define DYN_TYPE_RADGYR    22 
#define DYN_TYPE_DT_GAUSS  23 
#define DYN_TYPE_DT        24

#define DYN_TEMPERATURE 1
#define DYN_FORCE_CONS  2
#define DYN_ATOM_RADIUS 3
#define DYN_RADIUS_GYR  4
#define DYN_COORDS      5
#define DYN_ANGLE       6
#define DYN_TORSION     7
#define DYN_DIST        9
#define DYN_CHAIN       10 
#define DYN_SEG         11 
#define DYN_CHAIN_MODE  12 
#define DYN_SEG_MODE    13 
#define DYN_ATOMNAME    14
#define DYN_RESNAME     15
#define DYN_TENSOR      16
#define DYN_TEMPFACTOR  17

#define DYN_STORE     1
#define DYN_RESTORE   2

#define DYN_ALIGN     1
#define DYN_HOMOLOGY  2
#define DYN_SURFACE   3

#define DYN_PHIONLY   1
#define DYN_PSIONLY   2
#define DYN_PHIPSI    3

#define DYN_INFO_NULL       0
#define DYN_INFO_ATOMNAME   1
#define DYN_INFO_RESNAME    2
#define DYN_INFO_CHAINNAME  3
#define DYN_INFO_SEGNAME    4
#define DYN_INFO_LIMITS     5
#define DYN_INFO_EXTENT     6
#define DYN_INFO_DC         7
#define DYN_INFO_CS         8
#define DYN_INFO_PCS        9
#define DYN_INFO_NOE       10
#define DYN_INFO_J         11 
#define DYN_INFO_TORSION   12 
#define DYN_INFO_MW        13
#define DYN_INFO_RADGYR    14

#define DYN_MAP_SELECTED_ONLY 0 
#define DYN_MAP_ALL           1

#define DYN_TRANSFORM       1
#define DYN_TRANSLATE       2
#define DYN_ROTATE          3

#define DYN_BOLTZMAN      1.98719e-3 
#define DYN_PI            3.14159265358979323846264338

int dynInit(), dynInitTcl();
int dynCreate(), dynGet(), dynSet(), dynDestroy();
int dynRead(), dynWrite();

int dynSetInfo(), dynGetInfo(), dynGetMetric(), dynGetDist(), dynGetPhiPsi(), 
    dynGetPDBStr(), dynSetPhiPsi();

int dynGetAlpha(), dynGetKappa(), dynGetOmega(), dynSetOmega(), dynGetChi(), 
    dynSetChi(), dynGetNAAngles();

int dynXYZTrans(), dynGetCoord(), dynSetCoord();

int dynStore(), dynCopy(), dynSimulate(), dynSearch(), dynAlign(), dynSurface();
int dynSelect(), dynTransform(), dynEuler(), dynUpdateAC();
int dynEnergy(), dynDynamics();

int    dynSetTemperature();
double dynGetTemperature();
double dynGetKineticEnergy();
double dynTorsionAngle();

double dynSingleNOEDev();

int    dynSetZeroForces();
int    dynEnergyGrad2Force();

#define DYNRAND1 drand48
