/***/
/* hsymbol.h: definitions for symbol table use.
/***/

#define HSYM_NULL   0
#define HSYM_DOUBLE 1
#define HSYM_FLOAT  2
#define HSYM_INT    3
#define HSYM_STR    4
#define HSYM_PTR    5

#define HGET_D( HSYM ) HSYM->val.u.d
#define HGET_F( HSYM ) HSYM->val.u.f
#define HGET_I( HSYM ) HSYM->val.u.i
#define HGET_P( HSYM ) HSYM->val.u.s
#define HGET_S( HSYM ) HSYM->val.u.s

#define HSET_D( HSYM, VAL ) HSYM->val.u.d = VAL
#define HSET_F( HSYM, VAL ) HSYM->val.u.f = VAL
#define HSET_I( HSYM, VAL ) HSYM->val.u.i = VAL
#define HSET_P( HSYM, VAL ) HSYM->val.u.s = VAL
#define HSET_S( HSYM, VAL ) (void)strcpy( HSYM->val.u.s, VAL )

struct HTabInfo
   {
    int     hashLength;      /* Prime number hash length.              */
    int     symbolCount;     /* Number of valid symbols in this table. */
    int     maxDepth;        /* Max depth achieved by any symbol.      */
    struct  HSymInfo **base; /* First layer of hash table.             */
   };

struct HValInfo
   {
    int  type;        /* Symbol value type.                    */
    long int length;  /* Malloc symbol value length (strings). */

    union
       {
        double d;     /* Value for real symbols.               */
        float  f;     /* Value for real symbols.               */
        int    i;     /* Value for integer symbols.            */
        char   *s;    /* Value for string symbol (as pointer). */
       } u;
   };

struct HSymInfo
   {
    char   *name;            /* Symbol name.                               */
    int    depth;            /* Depth in symbol table.                     */
    int    id;               /* Unique ID.                                 */
    struct HValInfo  val;    /* Value associated with symbol.              */
    struct HSymInfo *next;   /* Pointer to next layer in table.            */ 
    struct HSymInfo *prev;   /* Pointer to prev layer in table.            */ 
    struct HTabInfo *parent; /* Pointer to symbol's parent table.          */ 
   };

typedef struct HTabInfo HTable;
typedef struct HSymInfo HSymbol;
typedef struct HValInfo HValue;

/***/
/* Prototypes for symbol table functions:
/***/

HTable *hSymTabAlloc();

int hValNull(), hSymNull(), hTabNull();
int hSymTabFree(), hSymFree(), hSymFreeR(), hSymDump();
int hSymExists(), hSymNameLoc();

HSymbol  *hSymAlloc(), *hFindOrAddSym(), *hFindSym(), *hFindSymID(), *hAddSym();
HSymbol  *hFindSymInt(), *hFindSymStr(), *hFindSymPtr();
char     *hNameAlloc(), *hSymVal2Str();

int hSetSymDlb(), hSetSymFlt(), hSetSymInt(), hSetSymStr(), hSetSymPtr();
int hGetSymDlb(), hGetSymFlt(), hGetSymInt(), hGetSymStr(), hGetSymPtr();
