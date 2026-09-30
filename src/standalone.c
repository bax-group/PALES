/*
 * Small, self-contained replacements for the NMRPipe/DYNAMO entry points
 * used by PALES 2.1.
 *
 * The copy of the DYNAMO sources distributed with this archive belongs to a
 * newer, much larger API and depends on Tcl plus many headers that are not
 * present.  PALES itself only needs command-line handling (cmndargs.c), the
 * allocator (memory.c), PDB I/O/selection, RDC table I/O/mapping, and the
 * five-parameter Saupe least-squares fit.  Those latter facilities live here
 * so that PALES can be built and used without an NMRPipe installation.
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prec.h"
#include "cmndargs.h"
#include "strmatch.h"
#include "hsymbol.h"
#include "pdb.h"
#include "dc.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define COMPAT_LINE_MAX 8192
#define COMPAT_COLUMN_MAX 64
#define COMPAT_NAME_MAX 32

/* cmndargs historically used an integer-returning strlen wrapper. */
int istrlen(const char *text)
{
    return text ? (int)strlen(text) : 0;
}

/* Used by the bundled case-insensitive regular-expression matcher. */
int raiseStr(char *text)
{
    char *p;

    if (!text)
        return 0;
    for (p = text; *p; ++p)
        *p = (char)toupper((unsigned char)*p);
    return 0;
}

static void copyTrim(char *out, size_t outSize, const char *in,
                     size_t start, size_t width)
{
    size_t length;
    size_t first;

    if (!out || outSize == 0)
        return;
    out[0] = '\0';
    if (!in)
        return;

    length = strlen(in);
    if (start >= length)
        return;
    if (start + width > length)
        width = length - start;
    first = 0;
    while (first < width && isspace((unsigned char)in[start + first]))
        ++first;
    while (width > first && isspace((unsigned char)in[start + width - 1]))
        --width;
    width -= first;
    if (width >= outSize)
        width = outSize - 1;
    memcpy(out, in + start + first, width);
    out[width] = '\0';
}

static void copyName(char *out, const char *in)
{
    size_t n;

    if (!out)
        return;
    if (!in) {
        out[0] = '\0';
        return;
    }
    n = strlen(in);
    if (n >= PDB_STRLEN)
        n = PDB_STRLEN - 1;
    memcpy(out, in, n);
    out[n] = '\0';
}

static int sameText(const char *a, const char *b)
{
    unsigned char ca;
    unsigned char cb;

    if (!a || !b)
        return 0;
    while (*a && *b) {
        ca = (unsigned char)toupper((unsigned char)*a++);
        cb = (unsigned char)toupper((unsigned char)*b++);
        if (ca != cb)
            return 0;
    }
    return *a == '\0' && *b == '\0';
}

static int isBlankName(const char *s)
{
    return !s || !*s || sameText(s, "NONE") || !strcmp(s, ".") ||
           !strcmp(s, "-");
}

int atomNull(struct AtomInfo *atom)
{
    int i;

    if (!atom)
        return 0;
    memset(atom, 0, sizeof(*atom));
    atom->shift = atom->dshift = atom->rcShift = 99999.0f;
    for (i = 0; i < DYN_NPREV; ++i)
        atom->time[i] = 0.0;
    return 0;
}

int pdbNull(struct PDBInfo *pdb)
{
    if (!pdb)
        return 0;
    memset(pdb, 0, sizeof(*pdb));
    pdb->firstAtomLine = -1;
    pdb->lastAtomLine = -1;
    return 0;
}

int pdbFree(struct PDBInfo *pdb)
{
    int i;

    if (!pdb)
        return 0;
    if (pdb->atomList) {
        for (i = 0; i < pdb->atomCount; ++i)
            free(pdb->atomList[i].auxInfo);
    }
    free(pdb->atomList);
    free(pdb->coordList);
    free(pdb->vdwList);
    free(pdb->resList);
    free(pdb->chainList);
    free(pdb->segList);
    free(pdb->bondMatrix);
    free(pdb->bondCount);
    free(pdb->neighborList);
    free(pdb->neighborCount);
    return pdbNull(pdb);
}

int pdbSelectNull(struct PDBSelectInfo *selection)
{
    if (!selection)
        return 0;
    memset(selection, 0, sizeof(*selection));
    selection->entry1 = PDB_BAD_RESID_VAL;
    selection->entryN = PDB_BAD_RESID_VAL;
    strcpy(selection->resNameList, "*");
    strcpy(selection->atomNameList, "*");
    strcpy(selection->chainNameList, "*");
    strcpy(selection->segNameList, "*");
    selection->resNameCount = selection->atomNameCount = 1;
    selection->chainNameCount = selection->segNameCount = 1;
    selection->resIDMode = STRMATCH_ALL;
    selection->resNameMode = STRMATCH_ALL;
    selection->atomNameMode = STRMATCH_EXACT;
    selection->chainNameMode = STRMATCH_ALL;
    selection->segNameMode = STRMATCH_ALL;
    selection->selectID = 1;
    selection->selectMode = PDB_SELECT_ALL;
    selection->selectClass = PDB_SELECT_ALL;
    selection->waterFlag = 1;
    selection->atomFlag = 1;
    return 0;
}

static int parsePDBAtom(const char *line, struct AtomInfo *atom, int model)
{
    char field[64];
    char record[8];
    size_t length;

    length = strlen(line);
    if (length < 54)
        return 0;
    copyTrim(record, sizeof(record), line, 0, 6);
    if (strcmp(record, "ATOM") && strcmp(record, "HETATM"))
        return 0;

    atomNull(atom);
    atom->hetFlag = !strcmp(record, "HETATM");
    atom->model = model;
    atom->selected = 1;

    copyTrim(field, sizeof(field), line, 6, 5);
    atom->atomID = (int)strtol(field, NULL, 10);
    copyTrim(atom->atomName, sizeof(atom->atomName), line, 12, 4);
    copyTrim(atom->locID, sizeof(atom->locID), line, 16, 1);
    copyTrim(atom->resName, sizeof(atom->resName), line, 17, 4);
    copyTrim(atom->chainID, sizeof(atom->chainID), line, 21, 1);
    copyTrim(field, sizeof(field), line, 22, 4);
    atom->resID = (int)strtol(field, NULL, 10);
    copyTrim(atom->icode, sizeof(atom->icode), line, 26, 1);
    copyTrim(field, sizeof(field), line, 30, 8);
    atom->x = (float)strtod(field, NULL);
    copyTrim(field, sizeof(field), line, 38, 8);
    atom->y = (float)strtod(field, NULL);
    copyTrim(field, sizeof(field), line, 46, 8);
    atom->z = (float)strtod(field, NULL);
    copyTrim(field, sizeof(field), line, 54, 6);
    atom->occupancy = *field ? (float)strtod(field, NULL) : 1.0f;
    copyTrim(field, sizeof(field), line, 60, 6);
    atom->tempFactor = *field ? (float)strtod(field, NULL) : 0.0f;
    copyTrim(atom->segID, sizeof(atom->segID), line, 72, 4);
    copyTrim(atom->element, sizeof(atom->element), line, 76, 2);
    copyTrim(atom->chargeStr, sizeof(atom->chargeStr), line, 78, 2);
    if (!atom->element[0]) {
        const char *p;
        p = atom->atomName;
        while (*p && isdigit((unsigned char)*p))
            ++p;
        if (*p) {
            atom->element[0] = (char)toupper((unsigned char)*p);
            atom->element[1] = '\0';
        }
    }
    return 1;
}

int readPDB(char *fileName, struct PDBInfo *pdb,
            struct PDBSelectInfo *selection)
{
    FILE *stream;
    struct AtomInfo atom;
    struct AtomInfo *newList;
    char line[COMPAT_LINE_MAX];
    int model;
    int capacity;
    int keep;

    if (!fileName || !*fileName || !pdb)
        return 1;
    stream = fopen(fileName, "r");
    if (!stream) {
        fprintf(stderr, "PALES: cannot open PDB file %s\n", fileName);
        return 1;
    }

    pdbFree(pdb);
    capacity = 0;
    model = 1;
    while (fgets(line, sizeof(line), stream)) {
        if (!strncmp(line, "MODEL", 5)) {
            model = (int)strtol(line + 5, NULL, 10);
            if (model < 1)
                model = 1;
            continue;
        }
        if (!parsePDBAtom(line, &atom, model))
            continue;
        keep = !atom.hetFlag;
        if (selection && selection->useRecordTypeFlag)
            keep = atom.hetFlag ? selection->hetatmFlag : selection->atomFlag;
        else if (selection && atom.hetFlag)
            keep = selection->hetatmFlag;
        if (!keep)
            continue;
        if (atom.locID[0] && atom.locID[0] != 'A')
            continue;
        if (pdb->atomCount == capacity) {
            capacity = capacity ? capacity * 2 : 1024;
            newList = (struct AtomInfo *)realloc(
                pdb->atomList, (size_t)capacity * sizeof(*newList));
            if (!newList) {
                fclose(stream);
                pdbFree(pdb);
                fprintf(stderr, "PALES: out of memory reading %s\n", fileName);
                return 1;
            }
            pdb->atomList = newList;
        }
        pdb->atomList[pdb->atomCount++] = atom;
        if (pdb->atomCount == 1) {
            pdb->minResID = pdb->maxResID = atom.resID;
            copyName(pdb->defChainName, atom.chainID);
            copyName(pdb->defSegName, atom.segID);
            pdb->firstAtomLine = 1;
        } else {
            if (atom.resID < pdb->minResID)
                pdb->minResID = atom.resID;
            if (atom.resID > pdb->maxResID)
                pdb->maxResID = atom.resID;
        }
    }
    fclose(stream);

    if (!pdb->atomCount) {
        fprintf(stderr, "PALES: no ATOM records found in %s\n", fileName);
        return 1;
    }
    pdb->lastAtomLine = pdb->atomCount;
    pdb->maxAtomCount = capacity;
    pdb->chainCount = 1;
    pdb->segCount = 1;
    pdb->modelCount = model;
    return 0;
}

int writePDB(char *fileName, struct PDBInfo *pdb,
             struct PDBSelectInfo *selection, char *remark)
{
    FILE *stream;
    struct AtomInfo *atom;
    int i;
    const char *record;
    const char *element;

    if (!fileName || !*fileName || !pdb)
        return 1;
    stream = fopen(fileName, "w");
    if (!stream) {
        fprintf(stderr, "PALES: cannot write PDB file %s\n", fileName);
        return 1;
    }
    if (remark && *remark)
        fprintf(stream, "REMARK %s\n", remark);
    for (i = 0; i < pdb->atomCount; ++i) {
        atom = pdb->atomList + i;
        if (selection && !atom->selected)
            continue;
        record = atom->hetFlag ? "HETATM" : "ATOM  ";
        element = atom->element[0] ? atom->element : "";
        fprintf(stream,
                "%-6s%5d %-4s%1s%-3s %1s%4d%1s   "
                "%8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n",
                record, atom->atomID, atom->atomName,
                atom->locID, atom->resName, atom->chainID,
                atom->resID, atom->icode, atom->x, atom->y, atom->z,
                atom->occupancy, atom->tempFactor, atom->segID,
                element, atom->chargeStr);
    }
    fputs("END\n", stream);
    fclose(stream);
    return 0;
}

int storePDBCoords(struct PDBInfo *pdb)
{
    int i;

    if (!pdb || !pdb->atomList)
        return 0;
    if (pdb->maxCoordCount < pdb->atomCount) {
        free(pdb->coordList);
        pdb->coordList = (struct Coord3D *)malloc(
            (size_t)pdb->atomCount * sizeof(*pdb->coordList));
        if (!pdb->coordList) {
            pdb->maxCoordCount = 0;
            return 1;
        }
        pdb->maxCoordCount = pdb->atomCount;
    }
    for (i = 0; i < pdb->atomCount; ++i) {
        pdb->coordList[i].x = pdb->atomList[i].x;
        pdb->coordList[i].y = pdb->atomList[i].y;
        pdb->coordList[i].z = pdb->atomList[i].z;
    }
    return 0;
}

int restorePDBCoords(struct PDBInfo *pdb)
{
    int i;

    if (!pdb || !pdb->atomList || !pdb->coordList ||
        pdb->maxCoordCount < pdb->atomCount)
        return 1;
    for (i = 0; i < pdb->atomCount; ++i) {
        pdb->atomList[i].x = pdb->coordList[i].x;
        pdb->atomList[i].y = pdb->coordList[i].y;
        pdb->atomList[i].z = pdb->coordList[i].z;
    }
    return 0;
}

int pdbSetSelection(struct PDBInfo *pdb, int selected)
{
    int i;

    if (!pdb)
        return 0;
    for (i = 0; i < pdb->atomCount; ++i)
        pdb->atomList[i].selected = selected;
    return pdb->atomCount;
}

static int matchGlob(const char *value, const char *pattern)
{
    int negate;
    int matched;
    int first;
    int last;

    while (*pattern) {
        if (*pattern == '*') {
            while (*pattern == '*')
                ++pattern;
            if (!*pattern)
                return 1;
            while (*value) {
                if (matchGlob(value, pattern))
                    return 1;
                ++value;
            }
            return 0;
        }
        if (!*value)
            return 0;
        if (*pattern == '?') {
            ++value;
            ++pattern;
            continue;
        }
        if (*pattern == '[') {
            ++pattern;
            negate = *pattern == '!' || *pattern == '^';
            if (negate)
                ++pattern;
            matched = 0;
            while (*pattern && *pattern != ']') {
                first = toupper((unsigned char)*pattern++);
                if (*pattern == '-' && pattern[1] && pattern[1] != ']') {
                    ++pattern;
                    last = toupper((unsigned char)*pattern++);
                    if (toupper((unsigned char)*value) >= first &&
                        toupper((unsigned char)*value) <= last)
                        matched = 1;
                } else if (toupper((unsigned char)*value) == first) {
                    matched = 1;
                }
            }
            if (*pattern == ']')
                ++pattern;
            if (matched == negate)
                return 0;
            ++value;
            continue;
        }
        if (toupper((unsigned char)*pattern) !=
            toupper((unsigned char)*value))
            return 0;
        ++pattern;
        ++value;
    }
    return *value == '\0';
}

static int matchOnePattern(const char *value, const char *pattern, int mode)
{
    if (!pattern || !*pattern || !strcmp(pattern, "*"))
        return 1;
    if (mode == STRMATCH_REGEX)
        return strMatch((char *)(value ? value : ""), (char *)pattern,
                        STRMATCH_REGEX, STRMATCH_CASE_NO);
    if (mode == STRMATCH_GLOB)
        return matchGlob(value ? value : "", pattern);
    return sameText(value ? value : "", pattern);
}

static int matchPatternList(const char *value, const char *list, int mode)
{
    char copy[NAMELEN + 1];
    char *token;

    if (mode == STRMATCH_ALL || !list || !*list || !strcmp(list, "*"))
        return 1;
    strncpy(copy, list, NAMELEN);
    copy[NAMELEN] = '\0';
    token = strtok(copy, " ,\t\r\n");
    while (token) {
        if (matchOnePattern(value, token, mode))
            return 1;
        token = strtok(NULL, " ,\t\r\n");
    }
    return 0;
}

int pdbSelect(struct PDBInfo *pdb, struct PDBSelectInfo *selection)
{
    struct AtomInfo *atom;
    int i;
    int count;
    int match;
    float dx;
    float dy;
    float dz;

    if (!pdb || !selection)
        return 0;
    if (selection->selectMode == PDB_SELECT_CLEAR)
        return pdbSetSelection(pdb, 0);
    if (selection->selectMode == PDB_SELECT_SET ||
        selection->selectMode == PDB_SELECT_ALL)
        return pdbSetSelection(pdb, selection->selectID);

    count = 0;
    for (i = 0; i < pdb->atomCount; ++i) {
        atom = pdb->atomList + i;
        match = 1;
        if (selection->selectMode == PDB_SELECT_SPHERE) {
            dx = atom->x - selection->x1;
            dy = atom->y - selection->y1;
            dz = atom->z - selection->z1;
            match = dx * dx + dy * dy + dz * dz <=
                    selection->rad * selection->rad;
        } else if (selection->selectMode == PDB_SELECT_REGION) {
            match = atom->x >= selection->x1 && atom->x <= selection->xn &&
                    atom->y >= selection->y1 && atom->y <= selection->yn &&
                    atom->z >= selection->z1 && atom->z <= selection->zn;
        } else {
            if (selection->resIDMode != STRMATCH_ALL)
                match = atom->resID >= selection->res1 &&
                        atom->resID <= selection->resN;
            if (match)
                match = matchPatternList(atom->atomName,
                    selection->atomNameList, selection->atomNameMode);
            if (match)
                match = matchPatternList(atom->resName,
                    selection->resNameList, selection->resNameMode);
            if (match)
                match = matchPatternList(atom->chainID,
                    selection->chainNameList, selection->chainNameMode);
            if (match)
                match = matchPatternList(atom->segID,
                    selection->segNameList, selection->segNameMode);
            if (match && selection->entry1 > PDB_BAD_RESID_VAL)
                match = i + 1 >= selection->entry1 &&
                        i + 1 <= selection->entryN;
            if (match && selection->useRecordTypeFlag)
                match = atom->hetFlag ? selection->hetatmFlag :
                                        selection->atomFlag;
            if (match && !selection->waterFlag && atom->hetFlag &&
                (sameText(atom->resName, "HOH") ||
                 sameText(atom->resName, "WAT")))
                match = 0;
        }
        if (selection->selectMode == PDB_SELECT_INV)
            match = !match;
        if (match) {
            atom->selected = selection->selectID;
            ++count;
        }
    }
    return count;
}

int pdbSelect1(struct PDBInfo *pdb, int res1, int resN,
               char *atomNames, int selectID, int matchMode)
{
    struct PDBSelectInfo selection;

    pdbSelectNull(&selection);
    selection.res1 = res1;
    selection.resN = resN;
    selection.resIDMode = STRMATCH_EXACT;
    selection.selectID = selectID;
    selection.atomNameMode = matchMode;
    selection.selectMode = PDB_SELECT_PARSE;
    strncpy(selection.atomNameList, atomNames, NAMELEN);
    selection.atomNameList[NAMELEN] = '\0';
    return pdbSelect(pdb, &selection);
}

int getPDBSelectArgs(int argc, char **argv,
                     struct PDBSelectInfo *selection)
{
    int n;
    int nameSet;
    int regionSet;

    if (!selection)
        return 0;
    n = 0;
    n += intArgD(argc, argv, "-res", &selection->res1);
    if (flagLoc(argc, argv, "-res")) {
        selection->resN = selection->res1;
        selection->resIDMode = STRMATCH_EXACT;
    }
    n += intArgD(argc, argv, "-r1", &selection->res1);
    n += intArgD(argc, argv, "-rN", &selection->resN);
    if (flagLoc(argc, argv, "-r1") || flagLoc(argc, argv, "-rN"))
        selection->resIDMode = STRMATCH_EXACT;
    n += intArgD(argc, argv, "-e1", &selection->entry1);
    n += intArgD(argc, argv, "-eN", &selection->entryN);
    n += intArgD(argc, argv, "-id", &selection->selectID);
    nameSet = strArgD(argc, argv, "-atomName", selection->atomNameList);
    nameSet += strArgD(argc, argv, "-atom", selection->atomNameList);
    if (nameSet) selection->atomNameMode = STRMATCH_EXACT;
    n += nameSet;
    nameSet = strArgD(argc, argv, "-resName", selection->resNameList);
    if (nameSet) selection->resNameMode = STRMATCH_EXACT;
    n += nameSet;
    nameSet = strArgD(argc, argv, "-chain", selection->chainNameList);
    if (nameSet) selection->chainNameMode = STRMATCH_EXACT;
    n += nameSet;
    nameSet = strArgD(argc, argv, "-seg", selection->segNameList);
    if (nameSet) selection->segNameMode = STRMATCH_EXACT;
    n += nameSet;
    if (flagLoc(argc, argv, "-exact"))
        selection->atomNameMode = selection->resNameMode =
        selection->chainNameMode = selection->segNameMode = STRMATCH_EXACT;
    if (flagLoc(argc, argv, "-glob"))
        selection->atomNameMode = selection->resNameMode =
        selection->chainNameMode = selection->segNameMode = STRMATCH_GLOB;
    if (flagLoc(argc, argv, "-regex"))
        selection->atomNameMode = selection->resNameMode =
        selection->chainNameMode = selection->segNameMode = STRMATCH_REGEX;
    if (flagLoc(argc, argv, "-aExact")) selection->atomNameMode = STRMATCH_EXACT;
    if (flagLoc(argc, argv, "-aGlob")) selection->atomNameMode = STRMATCH_GLOB;
    if (flagLoc(argc, argv, "-aRegex")) selection->atomNameMode = STRMATCH_REGEX;
    if (flagLoc(argc, argv, "-rExact")) selection->resNameMode = STRMATCH_EXACT;
    if (flagLoc(argc, argv, "-rGlob")) selection->resNameMode = STRMATCH_GLOB;
    if (flagLoc(argc, argv, "-rRegex")) selection->resNameMode = STRMATCH_REGEX;
    if (flagLoc(argc, argv, "-cExact")) selection->chainNameMode = STRMATCH_EXACT;
    if (flagLoc(argc, argv, "-cGlob")) selection->chainNameMode = STRMATCH_GLOB;
    if (flagLoc(argc, argv, "-cRegex")) selection->chainNameMode = STRMATCH_REGEX;
    if (flagLoc(argc, argv, "-sExact")) selection->segNameMode = STRMATCH_EXACT;
    if (flagLoc(argc, argv, "-sGlob")) selection->segNameMode = STRMATCH_GLOB;
    if (flagLoc(argc, argv, "-sRegex")) selection->segNameMode = STRMATCH_REGEX;
    if (flagLoc(argc, argv, "-useType")) selection->useRecordTypeFlag = 1;
    if (flagLoc(argc, argv, "-nouseType")) selection->useRecordTypeFlag = 0;
    if (flagLoc(argc, argv, "-atomrec")) selection->atomFlag = 1;
    if (flagLoc(argc, argv, "-noatomrec")) selection->atomFlag = 0;
    if (flagLoc(argc, argv, "-hetatm") || flagLoc(argc, argv, "-hetatmrec"))
        selection->hetatmFlag = 1;
    if (flagLoc(argc, argv, "-nohetatm") || flagLoc(argc, argv, "-nohetatmrec"))
        selection->hetatmFlag = 0;
    if (flagLoc(argc, argv, "-water")) selection->waterFlag = 1;
    if (flagLoc(argc, argv, "-nowater")) selection->waterFlag = 0;
    regionSet = 0;
    regionSet += fltArgD(argc, argv, "-x1", &selection->x1);
    regionSet += fltArgD(argc, argv, "-xn", &selection->xn);
    regionSet += fltArgD(argc, argv, "-y1", &selection->y1);
    regionSet += fltArgD(argc, argv, "-yn", &selection->yn);
    regionSet += fltArgD(argc, argv, "-z1", &selection->z1);
    regionSet += fltArgD(argc, argv, "-zn", &selection->zn);
    n += regionSet;
    if (fltArgD(argc, argv, "-near", &selection->rad)) {
        fltArgD(argc, argv, "-x", &selection->x1);
        fltArgD(argc, argv, "-y", &selection->y1);
        fltArgD(argc, argv, "-z", &selection->z1);
        selection->selectMode = PDB_SELECT_SPHERE;
        ++n;
    } else if (regionSet) {
        selection->selectMode = PDB_SELECT_REGION;
    }
    if (flagLoc(argc, argv, "-clear"))
        selection->selectMode = PDB_SELECT_CLEAR;
    else if (flagLoc(argc, argv, "-set"))
        selection->selectMode = PDB_SELECT_SET;
    else if (flagLoc(argc, argv, "-inv"))
        selection->selectMode = PDB_SELECT_INV;
    else if (n && selection->selectMode != PDB_SELECT_SPHERE &&
             selection->selectMode != PDB_SELECT_REGION)
        selection->selectMode = PDB_SELECT_PARSE;
    return n;
}

int dTensorNull(struct DTensor *tensor)
{
    if (!tensor)
        return 0;
    memset(tensor, 0, sizeof(*tensor));
    tensor->tensorMode = DC_SVD;
    tensor->eulerFlag = 1;
    tensor->daMin = tensor->drMin = -1.0e-2f;
    tensor->daMax = tensor->drMax = 1.0e-2f;
    tensor->psiMin = tensor->thetaMin = tensor->phiMin = -360.0f;
    tensor->psiMax = tensor->thetaMax = tensor->phiMax = 360.0f;
    tensor->gridStep = 5.0f;
    return 0;
}

int dTabNull(struct DTable *table)
{
    if (!table)
        return 0;
    memset(table, 0, sizeof(*table));
    dTensorNull(&table->dTensor);
    table->dMinRes = 0;
    table->dMaxRes = 0;
    table->nullChainName = 1;
    table->nullSegName = 1;
    return 0;
}

int dTabFree(struct DTable *table)
{
    if (!table)
        return 0;
    free(table->dObsList);
    free(table->dCalcList);
    free(table->workR);
    free(table->workI);
    return dTabNull(table);
}

int csTabNull(struct CSTable *table)
{
    if (!table)
        return 0;
    memset(table, 0, sizeof(*table));
    table->nullChainName = table->nullSegName = 1;
    table->nullPCSChainName = table->nullPCSSegName = 1;
    return 0;
}

int csTabFree(struct CSTable *table)
{
    if (!table)
        return 0;
    free(table->csList);
    free(table->pcsList);
    free(table->work);
    return csTabNull(table);
}

static void dInfoNull(struct DInfo *info)
{
    memset(info, 0, sizeof(*info));
    info->valid = 1;
    info->selected = 1;
    info->d = DC_UNDEF;
    info->dCalc = DC_UNDEF;
    info->di = 0.0f;
    info->w = 1.0f;
    info->sigma = 1.0f;
}

static int splitWords(char *line, char **words, int maxWords)
{
    int count;
    char *word;

    count = 0;
    word = strtok(line, " \t\r\n");
    while (word && count < maxWords) {
        words[count++] = word;
        word = strtok(NULL, " \t\r\n");
    }
    return count;
}

static int columnIndex(char names[][COMPAT_NAME_MAX], int count,
                       const char *name)
{
    int i;

    for (i = 0; i < count; ++i)
        if (sameText(names[i], name))
            return i;
    return -1;
}

static int anyColumn(char names[][COMPAT_NAME_MAX], int count,
                     const char *a, const char *b, const char *c)
{
    int i;

    i = a ? columnIndex(names, count, a) : -1;
    if (i < 0 && b)
        i = columnIndex(names, count, b);
    if (i < 0 && c)
        i = columnIndex(names, count, c);
    return i;
}

static const char *wordAt(char **words, int wordCount, int index)
{
    return index >= 0 && index < wordCount ? words[index] : "";
}

static int appendDInfo(struct DTable *table, struct DInfo *info)
{
    struct DInfo *observed;
    struct DInfo *calculated;
    int newCapacity;

    if (table->dCount == table->maxCount) {
        newCapacity = table->maxCount ? table->maxCount * 2 : 256;
        observed = (struct DInfo *)realloc(table->dObsList,
                    (size_t)newCapacity * sizeof(*observed));
        if (!observed)
            return 1;
        table->dObsList = observed;
        calculated = (struct DInfo *)realloc(table->dCalcList,
                    (size_t)newCapacity * sizeof(*calculated));
        if (!calculated)
            return 1;
        table->dCalcList = calculated;
        table->maxCount = newCapacity;
    }
    table->dObsList[table->dCount] = *info;
    table->dCalcList[table->dCount] = *info;
    table->dCalcList[table->dCount].d = DC_UNDEF;
    table->dCalcList[table->dCount].dCalc = DC_UNDEF;
    ++table->dCount;
    return 0;
}

int readDC(char *fileName, struct DTable *table)
{
    FILE *stream;
    char line[COMPAT_LINE_MAX];
    char work[COMPAT_LINE_MAX];
    char *words[COMPAT_COLUMN_MAX];
    char columns[COMPAT_COLUMN_MAX][COMPAT_NAME_MAX];
    struct DInfo info;
    int columnCount;
    int wordCount;
    int i;
    int index;
    struct DTensor savedTensor;

    if (!fileName || !*fileName || !table)
        return 1;
    stream = fopen(fileName, "r");
    if (!stream) {
        fprintf(stderr, "PALES: cannot open RDC table %s\n", fileName);
        return 1;
    }
    savedTensor = table->dTensor;
    dTabFree(table);
    table->dTensor = savedTensor;
    columnCount = 0;
    while (fgets(line, sizeof(line), stream)) {
        strncpy(work, line, sizeof(work) - 1);
        work[sizeof(work) - 1] = '\0';
        wordCount = splitWords(work, words, COMPAT_COLUMN_MAX);
        if (!wordCount || words[0][0] == '#')
            continue;
        if (sameText(words[0], "VARS")) {
            columnCount = wordCount - 1;
            for (i = 0; i < columnCount; ++i) {
                strncpy(columns[i], words[i + 1], COMPAT_NAME_MAX - 1);
                columns[i][COMPAT_NAME_MAX - 1] = '\0';
            }
            continue;
        }
        if (sameText(words[0], "FORMAT") || sameText(words[0], "DATA") ||
            sameText(words[0], "REMARK") || sameText(words[0], "TITLE"))
            continue;

        dInfoNull(&info);
        info.id = table->dCount + 1;
        info.count = 1;
        if (columnCount) {
            index = anyColumn(columns, columnCount, "RESID_I", "RESID1", "RESI_I");
            info.resIDI = (int)strtol(wordAt(words, wordCount, index), NULL, 10);
            index = anyColumn(columns, columnCount, "RESNAME_I", "RESNAME1", "RESN_I");
            copyName(info.resNameI, wordAt(words, wordCount, index));
            index = anyColumn(columns, columnCount, "ATOMNAME_I", "ATOMNAME1", "ATOM_I");
            copyName(info.atomNameI, wordAt(words, wordCount, index));
            index = anyColumn(columns, columnCount, "CHAINNAME_I", "CHAIN_I", "CHAIN1");
            copyName(info.chainNameI, wordAt(words, wordCount, index));
            index = anyColumn(columns, columnCount, "SEGNAME_I", "SEG_I", "SEG1");
            copyName(info.segNameI, wordAt(words, wordCount, index));

            index = anyColumn(columns, columnCount, "RESID_J", "RESID2", "RESI_J");
            info.resIDJ = (int)strtol(wordAt(words, wordCount, index), NULL, 10);
            index = anyColumn(columns, columnCount, "RESNAME_J", "RESNAME2", "RESN_J");
            copyName(info.resNameJ, wordAt(words, wordCount, index));
            index = anyColumn(columns, columnCount, "ATOMNAME_J", "ATOMNAME2", "ATOM_J");
            copyName(info.atomNameJ, wordAt(words, wordCount, index));
            index = anyColumn(columns, columnCount, "CHAINNAME_J", "CHAIN_J", "CHAIN2");
            copyName(info.chainNameJ, wordAt(words, wordCount, index));
            index = anyColumn(columns, columnCount, "SEGNAME_J", "SEG_J", "SEG2");
            copyName(info.segNameJ, wordAt(words, wordCount, index));

            index = anyColumn(columns, columnCount, "D_OBS", "D", "RDC");
            info.d = (float)strtod(wordAt(words, wordCount, index), NULL);
            index = anyColumn(columns, columnCount, "DI", "D_MAX", "DIPOLAR");
            if (index >= 0)
                info.di = (float)strtod(wordAt(words, wordCount, index), NULL);
            index = anyColumn(columns, columnCount, "DD", "SIGMA", "D_ERR");
            if (index >= 0)
                info.sigma = (float)strtod(wordAt(words, wordCount, index), NULL);
            index = anyColumn(columns, columnCount, "W", "WEIGHT", NULL);
            if (index >= 0)
                info.w = (float)strtod(wordAt(words, wordCount, index), NULL);
        } else {
            if (wordCount < 7)
                continue;
            info.resIDI = (int)strtol(words[0], NULL, 10);
            copyName(info.resNameI, words[1]);
            copyName(info.atomNameI, words[2]);
            info.resIDJ = (int)strtol(words[3], NULL, 10);
            copyName(info.resNameJ, words[4]);
            copyName(info.atomNameJ, words[5]);
            info.d = (float)strtod(words[6], NULL);
            if (wordCount > 7)
                info.sigma = (float)strtod(words[7], NULL);
            if (wordCount > 8)
                info.w = (float)strtod(words[8], NULL);
        }
        if (!info.atomNameI[0] || !info.atomNameJ[0] ||
            DC_UNDEF_TEST(info.d))
            info.valid = 0;
        if (info.sigma <= 0.0f)
            info.sigma = 1.0f;
        if (isBlankName(info.chainNameI)) info.chainNameI[0] = '\0';
        if (isBlankName(info.chainNameJ)) info.chainNameJ[0] = '\0';
        if (isBlankName(info.segNameI)) info.segNameI[0] = '\0';
        if (isBlankName(info.segNameJ)) info.segNameJ[0] = '\0';
        if (info.chainNameI[0] || info.chainNameJ[0]) table->nullChainName = 0;
        if (info.segNameI[0] || info.segNameJ[0]) table->nullSegName = 0;
        if (appendDInfo(table, &info)) {
            fclose(stream);
            dTabFree(table);
            fprintf(stderr, "PALES: out of memory reading %s\n", fileName);
            return 1;
        }
        if (table->dCount == 1) {
            table->dMinRes = info.resIDI < info.resIDJ ? info.resIDI : info.resIDJ;
            table->dMaxRes = info.resIDI > info.resIDJ ? info.resIDI : info.resIDJ;
        } else {
            if (info.resIDI < table->dMinRes) table->dMinRes = info.resIDI;
            if (info.resIDJ < table->dMinRes) table->dMinRes = info.resIDJ;
            if (info.resIDI > table->dMaxRes) table->dMaxRes = info.resIDI;
            if (info.resIDJ > table->dMaxRes) table->dMaxRes = info.resIDJ;
        }
    }
    fclose(stream);
    table->dcList = table->dObsList;
    table->dcCount = table->dCount;
    if (!table->dCount) {
        fprintf(stderr, "PALES: no RDC records found in %s\n", fileName);
        return 1;
    }
    return 0;
}

static int atomNameEqual(const char *pdbName, const char *tableName)
{
    if (sameText(pdbName, tableName))
        return 1;
    if ((sameText(pdbName, "H") && sameText(tableName, "HN")) ||
        (sameText(pdbName, "HN") && sameText(tableName, "H")))
        return 1;
    return 0;
}

static struct AtomInfo *findMappedAtom(struct PDBInfo *pdb, int resID,
                                       const char *resName,
                                       const char *atomName,
                                       const char *chainName,
                                       const char *segName)
{
    struct AtomInfo *atom;
    int i;

    if (!pdb)
        return NULL;
    for (i = 0; i < pdb->atomCount; ++i) {
        atom = pdb->atomList + i;
        if (atom->resID != resID || !atomNameEqual(atom->atomName, atomName))
            continue;
        if (resName && *resName && !sameText(atom->resName, resName))
            continue;
        if (chainName && *chainName && !sameText(atom->chainID, chainName))
            continue;
        if (segName && *segName && !sameText(atom->segID, segName))
            continue;
        return atom;
    }
    return NULL;
}

int mapDCInfo(struct PDBInfo *pdb, struct DTable *table)
{
    struct DInfo *obs;
    struct DInfo *calc;
    int i;
    int count;

    if (!pdb || !table)
        return 0;
    count = 0;
    for (i = 0; i < table->dCount; ++i) {
        obs = table->dObsList + i;
        calc = table->dCalcList + i;
        obs->atomI[0] = findMappedAtom(pdb, obs->resIDI, obs->resNameI,
            obs->atomNameI, obs->chainNameI, obs->segNameI);
        obs->atomJ[0] = findMappedAtom(pdb, obs->resIDJ, obs->resNameJ,
            obs->atomNameJ, obs->chainNameJ, obs->segNameJ);
        obs->atomICount = obs->atomI[0] ? 1 : 0;
        obs->atomJCount = obs->atomJ[0] ? 1 : 0;
        calc->atomI[0] = obs->atomI[0];
        calc->atomJ[0] = obs->atomJ[0];
        calc->atomICount = obs->atomICount;
        calc->atomJCount = obs->atomJCount;
        calc->valid = obs->valid && obs->atomI[0] && obs->atomJ[0];
        if (calc->valid && obs->selected) {
            ++count;
        }
    }
    table->validCount = count;
    return count;
}

int selectD(struct DInfo *list, int count, int first, int last, int mode)
{
    int i;
    int selected;
    int inI;
    int inJ;

    selected = 0;
    for (i = 0; i < count; ++i) {
        inI = list[i].resIDI >= first && list[i].resIDI <= last;
        inJ = list[i].resIDJ >= first && list[i].resIDJ <= last;
        if (mode == DC_EXCLUDE)
            list[i].selected = inI && inJ;
        else if (mode == DC_INCLUDE_ALL)
            list[i].selected = inI || inJ;
        else {
            int maxResidue;
            maxResidue = list[i].resIDI > list[i].resIDJ ?
                         list[i].resIDI : list[i].resIDJ;
            list[i].selected = maxResidue >= first && maxResidue <= last;
        }
        list[i].valid = list[i].valid && list[i].selected;
        if (list[i].valid)
            ++selected;
    }
    return selected;
}

static double gammaForAtom(const struct AtomInfo *atom)
{
    int c;

    if (!atom)
        return 0.0;
    c = atom->element[0] ? atom->element[0] : atom->atomName[0];
    if (isdigit((unsigned char)c) && atom->atomName[1])
        c = atom->atomName[1];
    switch (toupper((unsigned char)c)) {
    case 'H': return gH;
    case 'N': return gN;
    case 'C': return gC;
    case 'D': return gD;
    case 'F': return gF;
    case 'P': return gP;
    default: return 0.0;
    }
}

static float calculateDI(struct AtomInfo *a, struct AtomInfo *b)
{
    double dx;
    double dy;
    double dz;
    double distance;
    double ratio;
    double ga;
    double gb;

    if (!a || !b)
        return 0.0f;
    dx = a->x - b->x;
    dy = a->y - b->y;
    dz = a->z - b->z;
    distance = sqrt(dx * dx + dy * dy + dz * dz);
    if (distance < DC_MIN_DIST)
        return 0.0f;
    ga = gammaForAtom(a);
    gb = gammaForAtom(b);
    if (ga == 0.0 || gb == 0.0)
        return 0.0f;
    ratio = (ga * gb) / (gH * gN);
    return (float)(DIVAL_HN * ratio *
        pow(RAD_N_HN / distance, 3.0));
}

int buildDC(struct PDBInfo *pdb, struct DTable *table, int offset,
            int segmentFirst, float *weights, int weightCount)
{
    struct DInfo *obs;
    struct DInfo *calc;
    int i;
    int count;
    int weightIndex;

    if (!pdb || !table)
        return 0;
    count = 0;
    for (i = 0; i < table->dCount; ++i) {
        obs = table->dObsList + i;
        calc = table->dCalcList + i;
        calc->atomI[0] = findMappedAtom(pdb, obs->resIDI + offset,
            obs->resNameI, obs->atomNameI, obs->chainNameI, obs->segNameI);
        calc->atomJ[0] = findMappedAtom(pdb, obs->resIDJ + offset,
            obs->resNameJ, obs->atomNameJ, obs->chainNameJ, obs->segNameJ);
        calc->atomICount = calc->atomI[0] ? 1 : 0;
        calc->atomJCount = calc->atomJ[0] ? 1 : 0;
        calc->valid = obs->valid && obs->selected && calc->atomI[0] && calc->atomJ[0];
        if (!calc->valid)
            continue;
        calc->di = obs->di;
        if (fabs(calc->di) < 1.0e-12)
            calc->di = calculateDI(calc->atomI[0], calc->atomJ[0]);
        obs->di = calc->di;
        weightIndex = obs->resIDI - segmentFirst;
        if (weights && weightIndex >= 0 && weightIndex < weightCount)
            obs->w = calc->w = weights[weightIndex];
        ++count;
    }
    table->validCount = count;
    return count;
}

static int solve5(double a[5][5], double b[5], double result[5])
{
    int i;
    int j;
    int k;
    int pivot;
    double largest;
    double value;
    double factor;

    for (i = 0; i < 5; ++i) {
        pivot = i;
        largest = fabs(a[i][i]);
        for (j = i + 1; j < 5; ++j) {
            value = fabs(a[j][i]);
            if (value > largest) {
                largest = value;
                pivot = j;
            }
        }
        if (largest < 1.0e-20)
            return 1;
        if (pivot != i) {
            for (k = i; k < 5; ++k) {
                value = a[i][k];
                a[i][k] = a[pivot][k];
                a[pivot][k] = value;
            }
            value = b[i]; b[i] = b[pivot]; b[pivot] = value;
        }
        for (j = i + 1; j < 5; ++j) {
            factor = a[j][i] / a[i][i];
            for (k = i; k < 5; ++k)
                a[j][k] -= factor * a[i][k];
            b[j] -= factor * b[i];
        }
    }
    for (i = 4; i >= 0; --i) {
        value = b[i];
        for (j = i + 1; j < 5; ++j)
            value -= a[i][j] * result[j];
        result[i] = value / a[i][i];
    }
    return 0;
}

static int designRow(struct DInfo *info, double row[5])
{
    struct AtomInfo *a;
    struct AtomInfo *b;
    double x;
    double y;
    double z;
    double length;

    a = info->atomI[0];
    b = info->atomJ[0];
    if (!a || !b)
        return 1;
    x = a->x - b->x;
    y = a->y - b->y;
    z = a->z - b->z;
    length = sqrt(x * x + y * y + z * z);
    if (length < DC_MIN_DIST)
        return 1;
    x /= length; y /= length; z /= length;
    row[0] = z * z - 0.5 * x * x - 0.5 * y * y;
    row[1] = 0.5 * x * x - 0.5 * y * y;
    row[2] = 2.0 * x * y;
    row[3] = 2.0 * x * z;
    row[4] = 2.0 * y * z;
    return 0;
}

int fitD(struct PDBInfo *pdb, struct DTable *table, int verbose)
{
    double normal[5][5];
    double rhs[5];
    double solution[5];
    double row[5];
    double weight;
    double prediction;
    int i;
    int j;
    int k;
    int count;
    int fixed;

    (void)pdb;
    (void)verbose;
    if (!table)
        return 1;
    memset(normal, 0, sizeof(normal));
    memset(rhs, 0, sizeof(rhs));
    memset(solution, 0, sizeof(solution));
    count = 0;
    fixed = table->dTensor.tensorMode == DC_SAUPE ||
            table->dTensor.tensorMode == DC_STERIC ||
            table->dTensor.tensorMode == DC_ELECTROSTATIC ||
            table->dTensor.tensorMode == DC_FIXED_ALL;

    if (fixed) {
        for (i = 0; i < 5; ++i)
            solution[i] = table->dTensor.userSaupe[i] != 0.0f ?
                          table->dTensor.userSaupe[i] : table->dTensor.saupe[i];
    } else {
        for (i = 0; i < table->dCount; ++i) {
            if (!table->dObsList[i].valid || !table->dCalcList[i].valid ||
                designRow(table->dCalcList + i, row) ||
                fabs(table->dCalcList[i].di) < 1.0e-12)
                continue;
            weight = table->dObsList[i].w;
            if (table->dObsList[i].sigma > 0.0f)
                weight /= table->dObsList[i].sigma;
            weight *= weight;
            for (j = 0; j < 5; ++j) {
                rhs[j] += weight * row[j] * table->dCalcList[i].di *
                          table->dObsList[i].d;
                for (k = 0; k < 5; ++k)
                    normal[j][k] += weight * row[j] * row[k] *
                                    table->dCalcList[i].di *
                                    table->dCalcList[i].di;
            }
            ++count;
        }
        if (count < 5 || solve5(normal, rhs, solution)) {
            fprintf(stderr,
                    "PALES: at least five independent mapped RDCs are required for an SVD fit\n");
            return 1;
        }
    }

    for (i = 0; i < 5; ++i)
        table->dTensor.saupe[i] = (float)solution[i];
    for (i = 0; i < table->dCount; ++i) {
        if (!table->dCalcList[i].valid || designRow(table->dCalcList + i, row))
            continue;
        prediction = 0.0;
        for (j = 0; j < 5; ++j)
            prediction += row[j] * solution[j];
        prediction *= table->dCalcList[i].di;
        table->dCalcList[i].d = (float)prediction;
        table->dCalcList[i].dCalc = (float)prediction;
        table->dObsList[i].dCalc = (float)prediction;
    }
    return 0;
}

int getStatDC(struct DInfo *observed, struct DInfo *calculated, int count,
              float *rms, float *chi2, int *validCount)
{
    double sum;
    double chi;
    double difference;
    int n;
    int i;

    sum = 0.0;
    chi = 0.0;
    n = 0;
    for (i = 0; i < count; ++i) {
        if (!observed[i].valid || !calculated[i].valid ||
            DC_UNDEF_TEST(observed[i].d) || DC_UNDEF_TEST(calculated[i].d))
            continue;
        difference = observed[i].d - calculated[i].d;
        sum += difference * difference;
        if (observed[i].sigma > 0.0f)
            chi += difference * difference /
                   (observed[i].sigma * observed[i].sigma);
        ++n;
    }
    if (rms) *rms = n ? (float)sqrt(sum / n) : 0.0f;
    if (chi2) *chi2 = n ? (float)(chi / n) : 0.0f;
    if (validCount) *validCount = n;
    return n;
}

int rot2euler(float rotation[3][3], float *psiList,
              float *thetaList, float *phiList)
{
    double theta;
    double cosine;
    int i;

    for (i = 0; i < 8; ++i)
        psiList[i] = thetaList[i] = phiList[i] = 0.0f;
    theta = -asin(fmax(-1.0, fmin(1.0, rotation[2][0])));
    cosine = cos(theta);
    if (fabs(cosine) > 1.0e-8) {
        psiList[0] = (float)(atan2(rotation[2][1], rotation[2][2]) * 180.0 / M_PI);
        phiList[0] = (float)(atan2(rotation[1][0], rotation[0][0]) * 180.0 / M_PI);
    } else {
        psiList[0] = (float)(atan2(rotation[0][1], rotation[1][1]) * 180.0 / M_PI);
        phiList[0] = 0.0f;
    }
    thetaList[0] = (float)(theta * 180.0 / M_PI);
    return 1;
}

int euler2saupe(float *saupe, float psiDegrees, float thetaDegrees,
                float phiDegrees, float da, float dr)
{
    double psi;
    double theta;
    double phi;
    double sp;
    double cp;
    double st;
    double ct;
    double sf;
    double cf;
    double r[3][3];
    double diagonal[3];
    double matrix[3][3];
    int i;
    int j;
    int k;

    psi = psiDegrees * M_PI / 180.0;
    theta = thetaDegrees * M_PI / 180.0;
    phi = phiDegrees * M_PI / 180.0;
    sp = sin(psi); cp = cos(psi);
    st = sin(theta); ct = cos(theta);
    sf = sin(phi); cf = cos(phi);
    r[0][0] = ct * cf;
    r[0][1] = sp * st * cf - cp * sf;
    r[0][2] = cp * st * cf + sp * sf;
    r[1][0] = ct * sf;
    r[1][1] = sp * st * sf + cp * cf;
    r[1][2] = cp * st * sf - sp * cf;
    r[2][0] = -st;
    r[2][1] = sp * ct;
    r[2][2] = cp * ct;
    diagonal[0] = -da - 1.5 * dr;
    diagonal[1] = -da + 1.5 * dr;
    diagonal[2] = 2.0 * da;
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) {
            matrix[i][j] = 0.0;
            for (k = 0; k < 3; ++k)
                matrix[i][j] += r[i][k] * diagonal[k] * r[j][k];
        }
    saupe[0] = (float)matrix[2][2];
    saupe[1] = (float)(matrix[0][0] - matrix[1][1]);
    saupe[2] = (float)matrix[0][1];
    saupe[3] = (float)matrix[0][2];
    saupe[4] = (float)matrix[1][2];
    return 0;
}

int reportEuler(FILE *stream, struct DTensor *tensor)
{
    int i;

    if (!stream || !tensor || !tensor->eulerFlag || tensor->eulerCount < 1)
        return 0;
    for (i = 0; i < tensor->eulerCount && i < 8; ++i)
        fprintf(stream, "DATA EULER %8.3f %8.3f %8.3f\n",
                tensor->psiList[i], tensor->thetaList[i], tensor->phiList[i]);
    return tensor->eulerCount;
}

/* The PALES code calls these DYNAMO lifecycle/update hooks, but does not
 * consume any state from them. */
int dynInit(int keyMode, int ramaMode)
{
    (void)keyMode;
    (void)ramaMode;
    return 0;
}
int dynDone() { return 0; }
int surfNull() { return 0; }
int surfFree() { return 0; }
int dipolarNull() { return 0; }
int dipolarFree() { return 0; }
int fitAlignFree() { return 0; }
int updatePDBAngles(struct PDBInfo *pdb) { (void)pdb; return 0; }
int updatePDBShifts(struct PDBInfo *pdb) { (void)pdb; return 0; }
int updatePDBSurf() { return 0; }

