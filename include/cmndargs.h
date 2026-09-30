/***/
/* Definitions for use of command line argument extraction procedures.
/***/

#define CHECK_NONE  0
#define CHECK_FLAGS 1
#define CHECK_ALL   2
#define CHECK_DBG   3
#define CHECK_ENV   4

int  initArgs(), checkArgs();
char *getProgName();

char *safe_strcpy();

int isFlag(), flagLoc(), nextFlag(), nextFlagTerminator();
int intArg(), intArgD(), dblArgD(), fltArgD(), strArgD(), ptrArgD();
int logArgD(), notLogArgD();
int nmrIntArgD();
int envIsTrue();

float  fltArg();
int    intArgDN(), fltArgDN();

double dblArg();
char   *strArg(), *ptrArg(), *getNthArg();

int    charListArg();
int    *iListArg();
float  *rListArg();
char   **sListArg(), **strListArg(), **ptrListArg(), **copyArgv();

struct ZListInfo
 {
  char *name;
  int  zOff;
  int  aOff;
 };

int    isZFmt(), isAZFmt(), getFmtCount(), getZList(), getZListN(), getZCount();

int badFlag;
extern int badFlag;
