
/***/
/* Define variables for table locations.
/***/

#define DefineAtomLocVars \
   static int resIDLoc = -1, resNameLoc = -1, atomNameLoc = -1, \
              chainNameLoc = -1, segNameLoc = -1

#define DefineAtomILocVars \
   static int resIDILoc = -1, resNameILoc = -1, atomNameILoc = -1, \
              chainNameILoc = -1, segNameILoc = -1

#define DefineAtomJLocVars \
   static int resIDJLoc = -1, resNameJLoc = -1, atomNameJLoc = -1, \
              chainNameJLoc = -1, segNameJLoc = -1

#define DefineAtomKLocVars \
   static int resIDKLoc = -1, resNameKLoc = -1, atomNameKLoc = -1, \
              chainNameKLoc = -1, segNameKLoc = -1

#define DefineAtomLLocVars \
   static int resIDLLoc = -1, resNameLLoc = -1, atomNameLLoc = -1, \
              chainNameLLoc = -1, segNameLLoc = -1

#define DefineAtomListLocVars( N ) \
   static int resIDLoc[N], resNameLoc[N], atomNameLoc[N], \
              chainNameLoc[N], segNameLoc[N], iiNth

#define InitAtomListLocVars( N ) \
   for( iiNth = 0; iiNth < N; iiNth++ ) \
      { \
       resIDLoc[iiNth]      = -1; \
       resNameLoc[iiNth]    = -1; \
       atomNameLoc[iiNth]   = -1; \
       chainNameLoc[iiNth]  = -1; \
       segNameLoc[iiNth]    = -1; \
      } iiNth = 0

/***/
/* Extract locations of table variables.
/***/

#define GetAtomTabLoc( vI, vC ) \
  getAtomTabLocs( vI, vC, "",   \
                  &resIDLoc, &resNameLoc, &atomNameLoc, \
                  &chainNameLoc, &segNameLoc )

#define GetAtomTabILoc( vI, vC )  \
  getAtomTabLocs( vI, vC, "_I",   \
                  &resIDILoc, &resNameILoc, &atomNameILoc, \
                  &chainNameILoc, &segNameILoc )

#define GetAtomTabJLoc( vI, vC )  \
  getAtomTabLocs( vI, vC, "_J",   \
                  &resIDJLoc, &resNameJLoc, &atomNameJLoc, \
                  &chainNameJLoc, &segNameJLoc )

#define GetAtomTabKLoc( vI, vC )  \
  getAtomTabLocs( vI, vC, "_K",   \
                  &resIDKLoc, &resNameKLoc, &atomNameKLoc, \
                  &chainNameKLoc, &segNameKLoc )

#define GetAtomTabLLoc( vI, vC )  \
  getAtomTabLocs( vI, vC, "_L",   \
                  &resIDLLoc, &resNameLLoc, &atomNameLLoc, \
                  &chainNameLLoc, &segNameLLoc )

#define GetAtomTabNthLoc( vI, vC, vN, N )  \
  getAtomTabLocs( vI, vC, vN,              \
                  &resIDLoc[N], &resNameLoc[N], &atomNameLoc[N], \
                  &chainNameLoc[N], &segNameLoc[N] )

/***/
/* Extract atom specification from table entry.
/***/

#define GetAtomTabVal( vI, vC, vTab, vPtr, rPtr ) \
  getAtomTabVals( vI, vC, rPtr,                      \
                  resIDLoc, resNameLoc, atomNameLoc, \
                  chainNameLoc, segNameLoc,          \
                  &vPtr->resID, vPtr->resName, vPtr->atomName, \
                  vPtr->chainName, vPtr->segName );            \
  if (vPtr->chainName[0]) vTab->nullChainName = 0; \
  if (vPtr->segName[0])   vTab->nullSegName   = 0

#define GetAtomTabIVal( vI, vC, vTab, vPtr, rPtr ) \
  getAtomTabVals( vI, vC, rPtr,                         \
                  resIDILoc, resNameILoc, atomNameILoc, \
                  chainNameILoc, segNameILoc,           \
                  &vPtr->resIDI, vPtr->resNameI, vPtr->atomNameI, \
                  vPtr->chainNameI, vPtr->segNameI ); \
  if (vPtr->chainNameI[0]) vTab->nullChainName = 0;   \
  if (vPtr->segNameI[0])   vTab->nullSegName   = 0

#define GetAtomTabJVal( vI, vC, vTab, vPtr, rPtr ) \
  getAtomTabVals( vI, vC, rPtr,                         \
                  resIDJLoc, resNameJLoc, atomNameJLoc, \
                  chainNameJLoc, segNameJLoc,           \
                  &vPtr->resIDJ, vPtr->resNameJ, vPtr->atomNameJ, \
                  vPtr->chainNameJ, vPtr->segNameJ ); \
  if (vPtr->chainNameJ[0]) vTab->nullChainName = 0;   \
  if (vPtr->segNameJ[0])   vTab->nullSegName   = 0

#define GetAtomTabKVal( vI, vC, vTab, vPtr, rPtr ) \
  getAtomTabVals( vI, vC, rPtr,                         \
                  resIDKLoc, resNameKLoc, atomNameKLoc, \
                  chainNameKLoc, segNameKLoc,           \
                  &vPtr->resIDK, vPtr->resNameK, vPtr->atomNameK, \
                  vPtr->chainNameK, vPtr->segNameK ); \
  if (vPtr->chainNameK[0]) vTab->nullChainName = 0;   \
  if (vPtr->segNameK[0])   vTab->nullSegName   = 0

#define GetAtomTabLVal( vI, vC, vTab, vPtr, rPtr ) \
  getAtomTabVals( vI, vC, rPtr,                         \
                  resIDLLoc, resNameLLoc, atomNameLLoc, \
                  chainNameLLoc, segNameLLoc,           \
                  &vPtr->resIDL, vPtr->resNameL, vPtr->atomNameL, \
                  vPtr->chainNameL, vPtr->segNameL ); \
  if (vPtr->chainNameL[0]) vTab->nullChainName = 0;   \
  if (vPtr->segNameL[0])   vTab->nullSegName   = 0

#define GetAtomTabNthVal( vI, vC, vTab, vPtr, rPtr, N ) \
  getAtomTabVals( vI, vC, rPtr,                         \
                  resIDLoc[N], resNameLoc[N], atomNameLoc[N], \
                  chainNameLoc[N], segNameLoc[N],            \
                  &vPtr->resID[N], vPtr->resName[N], vPtr->atomName[N], \
                  vPtr->chainName[N], vPtr->segName[N] ); \
  if (vPtr->chainName[N][0]) vTab->nullChainName = 0;   \
  if (vPtr->segName[N][0])   vTab->nullSegName   = 0

/***/
/* Print FORMAT info for table output.
/***/

#define PrintAtomFmt( outUnit, vTab ) \
   FPR( outUnit, " %%4d %%4s %%4s" ); \
   if (!vTab->nullChainName) FPR( outUnit, " %%4s" ); \
   if (!vTab->nullSegName)   FPR( outUnit, " %%4s" )

#define PrintAtomIFmt( outUnit, vTab ) \
   FPR( outUnit, " %%4d %%4s %%4s" ); \
   if (!vTab->nullChainName) FPR( outUnit, " %%4s" ); \
   if (!vTab->nullSegName)   FPR( outUnit, " %%4s" )

#define PrintAtomJFmt( outUnit, vTab ) \
   FPR( outUnit, " %%4d %%4s %%4s" ); \
   if (!vTab->nullChainName) FPR( outUnit, " %%4s" ); \
   if (!vTab->nullSegName)   FPR( outUnit, " %%4s" )

#define PrintAtomKFmt( outUnit, vTab ) \
   FPR( outUnit, " %%4d %%4s %%4s" ); \
   if (!vTab->nullChainName) FPR( outUnit, " %%4s" ); \
   if (!vTab->nullSegName)   FPR( outUnit, " %%4s" )

#define PrintAtomLFmt( outUnit, vTab ) \
   FPR( outUnit, " %%4d %%4s %%4s" ); \
   if (!vTab->nullChainName) FPR( outUnit, " %%4s" ); \
   if (!vTab->nullSegName)   FPR( outUnit, " %%4s" )

/***/
/* Print VARS info for table output.
/***/

#define PrintAtomVars( outUnit, vTab ) \
   FPR( outUnit, " RESID RESNAME ATOMNAME" ); \
   if (!vTab->nullChainName) FPR( outUnit, " CHAINNAME" ); \
   if (!vTab->nullSegName)   FPR( outUnit, " SEGNAME" )

#define PrintAtomIVars( outUnit, vTab ) \
    FPR( outUnit, " RESID_I RESNAME_I ATOMNAME_I" ); \
    if (!vTab->nullChainName) FPR( outUnit, " CHAINNAME_I" ); \
    if (!vTab->nullSegName)   FPR( outUnit, " SEGNAME_I" )

#define PrintAtomJVars( outUnit, vTab ) \
    FPR( outUnit, " RESID_J RESNAME_J ATOMNAME_J" ); \
    if (!vTab->nullChainName) FPR( outUnit, " CHAINNAME_J" ); \
    if (!vTab->nullSegName)   FPR( outUnit, " SEGNAME_J" )

#define PrintAtomKVars( outUnit, vTab ) \
    FPR( outUnit, " RESID_K RESNAME_K ATOMNAME_K" ); \
    if (!vTab->nullChainName) FPR( outUnit, " CHAINNAME_K" ); \
    if (!vTab->nullSegName)   FPR( outUnit, " SEGNAME_K" )

#define PrintAtomLVars( outUnit, vTab ) \
    FPR( outUnit, " RESID_L RESNAME_L ATOMNAME_L" ); \
    if (!vTab->nullChainName) FPR( outUnit, " CHAINNAME_L" ); \
    if (!vTab->nullSegName)   FPR( outUnit, " SEGNAME_L" )

#define PrintAtomNthVars( outUnit, vTab, N ) \
    FPR( outUnit, " RESID_%c RESNAME_%c ATOMNAME_%c", 'A' + N ); \
    if (!vTab->nullChainName) FPR( outUnit, " CHAINNAME_%c", 'A' + N ); \
    if (!vTab->nullSegName)   FPR( outUnit, " SEGNAME_%c",   'A' + N )

/***/
/* Print atom info for table output.
/***/

#define PrintAtomVals( outUnit, vPtr, vTab ) \
 FPR( outUnit, " %4d %4s %4s", vPtr->resID, vPtr->resName, vPtr->atomName ); \
 if (!vTab->nullChainName) \
    FPR( outUnit, " %4s", vPtr->chainName[0] ? vPtr->chainName : "None" );   \
 if (!vTab->nullSegName)   \
    FPR( outUnit, " %4s", vPtr->segName[0]   ? vPtr->segName   : "None" )

#define PrintAtomIVals( outUnit, vPtr, vTab ) \
 FPR(outUnit, " %4d %4s %4s", vPtr->resIDI, vPtr->resNameI, vPtr->atomNameI); \
 if (!vTab->nullChainName) \
    FPR( outUnit, " %4s", vPtr->chainNameI[0] ? vPtr->chainNameI : "None" );  \
 if (!vTab->nullSegName)   \
    FPR( outUnit, " %4s", vPtr->segNameI[0]   ? vPtr->segNameI   : "None" )

#define PrintAtomJVals( outUnit, vPtr, vTab ) \
 FPR(outUnit, " %4d %4s %4s", vPtr->resIDJ, vPtr->resNameJ, vPtr->atomNameJ); \
 if (!vTab->nullChainName) \
    FPR( outUnit, " %4s", vPtr->chainNameJ[0] ? vPtr->chainNameJ : "None" );  \
 if (!vTab->nullSegName)   \
    FPR( outUnit, " %4s", vPtr->segNameJ[0]   ? vPtr->segNameJ   : "None" )

#define PrintAtomKVals( outUnit, vPtr, vTab ) \
 FPR(outUnit, " %4d %4s %4s", vPtr->resIDK, vPtr->resNameK, vPtr->atomNameK); \
 if (!vTab->nullChainName) \
    FPR( outUnit, " %4s", vPtr->chainNameK[0] ? vPtr->chainNameK : "None" );  \
 if (!vTab->nullSegName)   \
    FPR( outUnit, " %4s", vPtr->segNameK[0]   ? vPtr->segNameK   : "None" )

#define PrintAtomLVals( outUnit, vPtr, vTab ) \
 FPR(outUnit, " %4d %4s %4s", vPtr->resIDL, vPtr->resNameL, vPtr->atomNameL); \
 if (!vTab->nullChainName) \
    FPR( outUnit, " %4s", vPtr->chainNameL[0] ? vPtr->chainNameL : "None" );  \
 if (!vTab->nullSegName)   \
    FPR( outUnit, " %4s", vPtr->segNameL[0]   ? vPtr->segNameL   : "None" )

#define PrintAtomNthVals( outUnit, vPtr, vTab, N ) \
FPR(outUnit," %4d %4s %4s",vPtr->resID[N],vPtr->resName[N],vPtr->atomName[N]); \
if (!vTab->nullChainName) \
   FPR( outUnit, " %4s", vPtr->chainName[N][0] ? vPtr->chainName[N]:"None" );  \
if (!vTab->nullSegName)   \
   FPR( outUnit, " %4s", vPtr->segName[N][0]   ? vPtr->segName[N]  :"None" )

#define PrintAtomNullVals( outUnit, vPtr, vTab ) \
 FPR(outUnit, " %4d %4s %4s", 0, "None", "None"); \
 if (!vTab->nullChainName) \
    FPR( outUnit, " %4s", "None" );  \
 if (!vTab->nullSegName)   \
    FPR( outUnit, " %4s", "None" )

/***/
/* Adjustments for unspecified chain and segment:
/***/

#define DefChain( pdbInfo, vTab, name ) \
   vTab->nullChainName ? pdbInfo->defChainName : name

#define DefSeg( pdbInfo, vTab, name ) \
   vTab->nullSegName   ? pdbInfo->defSegName   : name

/***/
/* Print Atom Spec for Warning Messages:
/***/

#define ShowAtomInfo( fUnit, vPtr ) \
 FPR( fUnit, \
      "   ResID: %4d ResName: %4s AtomName: %4s Seg: %4s Chain: %4s\n", \
      vPtr->resID, vPtr->resName, vPtr->atomName, \
      vPtr->segID, vPtr->chainID )

#define ShowAtom( fUnit, vPtr ) \
 FPR( fUnit, \
      "   ResID: %4d ResName: %4s AtomName: %4s Seg: %4s Chain: %4s\n", \
      vPtr->resID, vPtr->resName, vPtr->atomName, \
      vPtr->segName, vPtr->chainName )

#define ShowAtomI( fUnit, vPtr ) \
 FPR( fUnit, \
      "   ResIDI: %4d ResNameI: %4s AtomNameI: %4s SegI: %4s ChainI: %4s\n", \
      vPtr->resIDI, vPtr->resNameI, vPtr->atomNameI, \
      vPtr->segNameI, vPtr->chainNameI )

#define ShowAtomJ( fUnit, vPtr ) \
 FPR( fUnit, \
      "   ResIDJ: %4d ResNameJ: %4s AtomNameJ: %4s SegJ: %4s ChainJ: %4s\n", \
      vPtr->resIDJ, vPtr->resNameJ, vPtr->atomNameJ, \
      vPtr->segNameJ, vPtr->chainNameJ )

#define ShowAtomK( fUnit, vPtr ) \
 FPR( fUnit, \
      "   ResIDK: %4d ResNameK: %4s AtomNameK: %4s SegK: %4s ChainK: %4s\n", \
      vPtr->resIDK, vPtr->resNameK, vPtr->atomNameK, \
      vPtr->segNameK, vPtr->chainNameK )

#define ShowAtomL( fUnit, vPtr ) \
 FPR( fUnit, \
      "   ResIDL: %4d ResNameL: %4s AtomNameL: %4s SegL: %4s ChainL: %4s\n", \
      vPtr->resIDL, vPtr->resNameL, vPtr->atomNameL, \
      vPtr->segNameL, vPtr->chainNameL )

#define ShowAtomNth( fUnit, vPtr, N ) \
FPR( fUnit, \
   "   ResID%c: %4d ResName%c: %4s AtomName%c: %4s Seg%c: %4s Chain%c: %4s\n", \
   'A'+N, vPtr->resID[N], 'A'+N, vPtr->resName[N], 'A'+N, vPtr->atomName[N],   \
   'A'+N, vPtr->segName[N], 'A'+N, vPtr->chainName[N] )

/***/
/* Parsing DATA SEQUENCE statements:
/***/

struct SeqInfo
   {
    int   firstResID;
    int   resCount;
    char  **resNameList;
   };

