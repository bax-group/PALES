/***/
/* NSUBEXP must be at least 10, and no greater than 117 or the parser
/* will not work properly.
/***/

#define STRMATCH_CASE_YES 1
#define STRMATCH_CASE_NO  0

#define STRMATCH_ALL      0
#define STRMATCH_EXACT    1
#define STRMATCH_REGEX    2
#define STRMATCH_GLOB     3

#define NSUBEXP  20

typedef struct Regexp
   {
    char *startp[NSUBEXP];
    char *endp[NSUBEXP];
    char regstart;
    char reganch;
    char *regmust;
    int regmlen;
    char program[1];
   } Regexp;

int    strRegExec();
int    strRegExec2();
void   strRegSub();
void   strRegError();

Regexp *strRegComp();
char   *strGetRegError();

int strMatch(), strMatchGlob(), strMatchRegex(), strbegin();
int strMaxIdentity();
int atomMatch();

