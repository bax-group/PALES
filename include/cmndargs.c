/******************************************************************************/
/*                                                                            */
/*                   ---- NIH NMR Software System ----                        */
/*                        Copyright 1992 and 1993                             */
/*                             Frank Delaglio                                 */
/*                   NIH Laboratory of Chemical Physics                       */
/*                                                                            */
/*               This software is not for distribution without                */
/*                  the written permission of the author.                     */
/*                                                                            */
/******************************************************************************/

/***/
/* cmndargs: tools for extracting arguments and lists from
/*           command line arguments.
/*
/* nullArgs:  null-initializes extraction and checking of arguments.
/* initArgs:  initializes extraction and checking of arguments.
/* checkArgs: performs final checking of arguments.
/*
/* isFlag:   returns true if argument is a flag.
/*
/* flagLoc:            returns location of flag in argument list.
/* nextFlag:           returns location of next flag in the argument list.
/* nextFlagTerminator: returns location of next -- in the argument list.
/*
/* intArg:   returns integer value associated with flag.
/* fltArg:   returns float value assosciated with flag.
/* strArg:   returns string pointer (buffer copy) associated with flag.
/* ptrArg:   returns string pointer (original) associated with flag.
/*
/* intArgD:  sets integer value if flag is present.
/* fltArgD:  sets float value if flag is present.
/* strArgD:  sets (copies) string value if flag is present.
/* ptrArgD:  sets pointer value if flag is present.
/*
/* intArgDN: returns a list of at most N integers in the existing array.
/* fltArgDN: returns a list of at most N floats in the existing array.
/*
/* iListArg:    allocates and returns integer list associated with flag. 
/* rListArg:    allocates and returns real list associated with flag. 
/* strListArg:  allocates and returns a string list associated with flag.
/* ptrListArg:  returns pointer value and list length associated with flag.
/* charListArg: returns list associated with a flag in a single string.
/*
/* envIsTrue:   return 1 if environment variable is defined as 1 Yes or True.
/***/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>

#include "cmndargs.h"
#include "memory.h"
#include "token.h"
#include "prec.h"
#include "atof.h"
#include "lintfix.h"

#define MAXARGS NAMELEN

static char *currentProgName = (char *)NULL;
static char argMap[MAXARGS+1];

static int  currentCheckMode = CHECK_NONE;
static int  currentArgc      = 0;
static int  currentChecksum  = 0;

static int  isBlank();

static char DEF_PROGRAM_NAME[] = "program";

#define FPR (void)fprintf

/***/
/* nullArgs: null-initialize variables for argument extraction.
/***/

int nullArgs()
{
    int loc;

    for( loc = 0; loc < MAXARGS; loc++ ) argMap[loc] = '0';

    argMap[MAXARGS]  = '\0';
    currentProgName  = (char *)NULL;
    currentCheckMode = CHECK_NONE;
    currentArgc      = 0;
    currentChecksum  = 0;

    return( 0 );
}

/***/
/* getArgvChecksum: finds a checksum of characters in argv.
/***/

int getArgvChecksum( argc, argv )

   int  argc;
   char **argv;
{
    int  i, sum;
    char *sPtr;

    sum = 0;

    for( i = 0; i < argc; i++ )
       {
        for( sPtr = argv[i]; sPtr && *sPtr; sPtr++ )
           {
            sum += (unsigned int) *sPtr;
           }
       }

    return( sum );
}

/***/
/* initArgs: initialize argument extraction and checking. 
/***/

int initArgs( argc, argv, checkMode )

   int  argc, checkMode;
   char **argv;
{
    char *sPtr;
    int  loc;

/***/
/* Initialize variables which record which args have been used:
/***/

    (void) nullArgs();

    argMap[0]        = '1';
    currentProgName  = argv[0];
    currentArgc      = argc;
    currentCheckMode = checkMode;

/***/
/* Test environment variable if needed to see what checking to perform:
/***/

    if (checkMode == CHECK_ENV)
       {
        if ((sPtr = getenv( "NMR_CHECKARGS" )))
           {
            if (!strcasecmp( sPtr, "ALL" ))
               currentCheckMode = CHECK_ALL;
            else if (!strcasecmp( sPtr, "FLAGS" ))
               currentCheckMode = CHECK_FLAGS;
            else if (!strcasecmp( sPtr, "DEBUG" ))
               currentCheckMode = CHECK_DBG;
            else
               currentCheckMode = CHECK_ALL;
           }
        else
           {
            currentCheckMode = CHECK_NONE;
           }
       }

/***/
/* If any checking will be used, record checksum of argv for later testing:
/***/

    if (currentCheckMode == CHECK_NONE) return( 0 );

    if (flagLoc( argc, argv, "-help" ))
       {
        FPR( stderr,
             "Program %s, argument checking enabled.\n",
             getProgName() );
       }

    if (argc > MAXARGS)
       {
        FPR( stderr,
             "Warning from %s: not all arguments will be checked.\n",
             getProgName() );
       }

    currentChecksum = getArgvChecksum( argc, argv );

    return( 0 );
}

/***/
/* checkArgs: check for unused arguments in command line.
/***/

int checkArgs( argc, argv )

   int  argc;
   char **argv;
{
    int badArg, loc1, loc2, checksum, count;

/***/
/* Stop if no checking is required:
/***/

    count = 0;

    if ((currentCheckMode == CHECK_NONE) || (argc < 2))
       {
        (void) nullArgs();
        return( count );
       }

/***/
/* Special debug mode: print map of used arguments:
/***/

    checksum = getArgvChecksum( argc, argv );

    if (currentCheckMode == CHECK_DBG)
       {
        FPR( stderr,
             "Debug Info from %s. Original Checksum: %d Final Checksum: %d.\n",
             getProgName(), currentChecksum, checksum );

        currentCheckMode = CHECK_ALL;

        for( loc1 = 0; (loc1 < argc) && (loc1 < MAXARGS); loc1++ )
           {
            FPR( stderr,
                 "%3d. %c %3d Characters: '%s'\n",
                 loc1, argMap[loc1],
                 argv[loc1] ? istrlen( argv[loc1] ) : (int)-1, 
                 argv[loc1] ? argv[loc1] : "(None)" );
           }
       }

/***/
/* Compare checksum or argv to the initial value.
/* Reduce argc to avoid overflow of argument mask.
/* Reset entries for blank argv elements, which are always unused. 
/***/

    if ((checksum != currentChecksum) || (argc != currentArgc))
       {
        FPR( stderr,
             "Warning from %s: argument checking may not be valid.\n",
             getProgName() );
       }

    argc = argc < MAXARGS ? argc : MAXARGS;

    for( loc1 = 1; loc1 < argc; loc1++ )
       {
        if (isBlank( argv[loc1] )) argMap[loc1] = '1';
       }

/***/
/* Find and display ranges of arguments which were unused.
/***/

    for( loc1 = 1; loc1 < argc; loc1++ )
       {
        if (currentCheckMode == CHECK_ALL)
           badArg = (argMap[loc1] == '0');
        else
           badArg = (argMap[loc1] == '0') && isFlag( argv[loc1] );

        if (badArg)
           {
            FPR( stderr, "Warning from %s:\n", getProgName() );
       
            if (currentCheckMode == CHECK_FLAGS)
               {
                loc2 = loc1;
               }
            else
               {     
                for( loc2 = loc1; loc2 < argc - 1; loc2++ )
                   {
                    if (argMap[loc2+1] != '0') break; 
                   }
               }

            if (loc1 != loc2)
               {
                FPR( stderr,
                     "  Arguments %d to %d may be unknown or unused:\n",
                     loc1, loc2 );
               }
            else
               {
                FPR( stderr,
                     "  Argument %d may be unknown or unused:\n",
                     loc1 );
               }

            FPR( stderr, "  '" );

            while( loc1 <= loc2 )
               {
                FPR( stderr, " %s", argv[loc1] );
                loc1++;
               }

            FPR( stderr, " '\n" );
            count++;
           }
       }

    (void) nullArgs();

    return( count );
}

/***/
/* getProgName: extract program name recorded by initArgs.
/***/

char *getProgName()
{
    if (!currentProgName)
       {
        return( DEF_PROGRAM_NAME );
       }

    return( currentProgName );
}

/***/
/* isNumber: return 1 if string is a valid number string.
/***/

int isNumber( string )

   char *string;
{
    float r;

    if (!string) return( 0 );

    if (1 == sscanf( string, "%f", &r ))
       return( 1 );
    else
       return( 0 );
}

/***/
/* isFlag: return true if argument is a flag.
/*         Currently, must begin with "-" character.
/***/

int isFlag( arg )

   char *arg;
{
   if (!arg) return( 0 );

   if (arg[0] == '-')
      {
       if (isdigit( arg[1] ) || arg[1] == '.')
          return( 0 );
       else
          return( 1 );
      }

    return( 0 );
}

/***/
/* flagLoc: return 0-origin location of flag in argv list.
/***/

int flagLoc( argc, argv, flag )

   int argc;
   char **argv, *flag;
{
    int loc;

    if (argc < 0 || !argv) 
       {
        return( 0 );
        badFlag = 1;
       }

    badFlag = 0;

    for( loc = 0; loc < argc; loc++ )
       {
        if (!strcmp( argv[loc], flag ))
           {
            if (loc < MAXARGS) argMap[loc] = '1';
            return( loc );
           }
       }

    badFlag = 1;

    return( 0 );
}

/***/
/* nextFlag: return location of next flag, if any, in argv list.
/***/

int nextFlag( argc, argv, prevFlag )

   int  argc, prevFlag;
   char **argv;
{
    int loc;

    badFlag = 0;

    for( loc = prevFlag + 1; loc < argc; loc++ )
       {
        if (isFlag( argv[loc] )) return( loc );
       }

    return( argc );
}

int nextFlagTerminator( argc, argv, prevFlag )

   int  argc, prevFlag;
   char **argv;
{
    int loc;

    badFlag = 0;

    for( loc = prevFlag + 1; loc < argc; loc++ )
       {
        if (!strcmp( argv[loc], "--" )) return( loc );
       }

    return( argc );
}

/***/
/* intArg: extracts an integer flag argument.
/***/

int intArg( argc, argv, flag )

   int  argc;
   char **argv, *flag;
{
   int loc;

   loc = 1 + flagLoc( argc, argv, flag );

   if (badFlag || loc < 1 || loc > argc - 1)
      {
       badFlag = 1;
       return( 0 );
      }
   else
      {
       badFlag = 0;

       if (loc < MAXARGS) argMap[loc] = '1';

       return( atoi( argv[loc] ));
      }
}

/***/
/* fltArg: extracts a floating point flag argument.
/***/

float fltArg( argc, argv, flag )

   int  argc;
   char **argv, *flag;
{
   int   loc;

   loc = 1 + flagLoc( argc, argv, flag );

   if (badFlag || loc < 1 || loc > argc - 1)
      {
       badFlag = 1;
       return( 0.0 );
      }
   else
      {
       badFlag = 0;

       if (loc < MAXARGS) argMap[loc] = '1';

       return( (float) ATOF( argv[loc] ));
      }
}

/***/
/* dblArg: extracts a floating point flag argument.
/***/

double dblArg( argc, argv, flag )

   int  argc;
   char **argv, *flag;
{
   int   loc;

   loc = 1 + flagLoc( argc, argv, flag );

   if (badFlag || loc < 1 || loc > argc - 1)
      {
       badFlag = 1;
       return( (double)0.0 );
      }
   else
      {
       badFlag = 0;

       if (loc < MAXARGS) argMap[loc] = '1';

       return( ATOF( argv[loc] ));
      }
}

/***/
/* strArg: extracts pointer to a copy of a string flag argument.
/*         String contents should be copied by caller.
/***/

char *strArg( argc, argv, flag )

   int  argc;
   char **argv, *flag;
{
   int    loc;
   static char string[NAMELEN+1];

   loc = 1 + flagLoc( argc, argv, flag );

   if (badFlag || loc < 1 || loc > argc - 1)
      {
       badFlag = 1;
       return( (char *) NULL );
      }
   else
      {
       badFlag = 0;

       if (loc < MAXARGS) argMap[loc] = '1';

       (void) strcpy( string, argv[loc] );
       return( string );
      }
}

/***/
/* getNthArg: returns argument n of argv list.
/***/

char *getNthArg( argc, argv, n )

   int  argc, n;
   char **argv;
{
    if (n < 0 || n > argc - 1) return( (char *)NULL );

    if (n < MAXARGS) argMap[n] = '1';

    return( argv[n] );
}

/***/
/* ptrArg: extracts pointer to a string flag argument.
/***/

char *ptrArg( argc, argv, flag )

   int  argc;
   char **argv, *flag;
{
   int    loc;

   loc = 1 + flagLoc( argc, argv, flag );

   if (badFlag || loc < 1 || loc > argc - 1)
      {
       badFlag = 1;
       return( (char *) NULL );
      }
   else
      {
       badFlag = 0;

       if (loc < MAXARGS) argMap[loc] = '1';

       return( argv[loc] );
      }
}

/***/
/* intArgDN: returns a list of at most n integers in the existing array;
/*           accomodates multiple items within each argument in argv.
/***/

int intArgDN( argc, argv, flag, iList, nMax )

   int  argc, *iList, nMax;
   char **argv, *flag;
{
   int i, n, loc, loc1, loc2, count;

   if (!iList || nMax < 1) return( -1 );

   loc1 = flagLoc( argc, argv, flag );
   if (badFlag) return( -1 );

   loc2 = nextFlag( argc, argv, loc1 );
   if (badFlag) return( -1 );

   count = 0;

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       if (loc < MAXARGS) argMap[loc] = '1';

       n = cntToken( argv[loc] );

       for( i = 1; i <= n; i++ )
          {
           if (count >= nMax) break;
           iList[count++] = getTokenI( i, argv[loc] );
          }
      }

   return( count );
}

/***/
/* fltArgDN: returns a list of at most n integers in the existing array;
/*           accomodates multiple items within each argument in argv.
/***/

int fltArgDN( argc, argv, flag, rList, nMax )

   float *rList;
   int   argc, nMax;
   char  **argv, *flag;
{
   int i, n, loc, loc1, loc2, count;

   if (!rList || nMax < 1) return( -1 );

   loc1 = flagLoc( argc, argv, flag );
   if (badFlag) return( -1 );

   loc2 = nextFlag( argc, argv, loc1 );
   if (badFlag) return( -1 );

   count = 0;

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       if (loc < MAXARGS) argMap[loc] = '1';

       n = cntToken( argv[loc] );

       for( i = 1; i <= n; i++ )
          {
           if (count >= nMax) break;
           rList[count++] = getTokenR( i, argv[loc] );
          }
      }

   return( count );
}

/***/
/* iListArg: allocates and returns an integer argument flag list;
/*           accomodates multiple items within each argument in argv.
/***/

int *iListArg( argc, argv, flag, length )

   int  argc, *length;
   char **argv, *flag;
{
   int i, n, loc, loc1, loc2, *iArray, *iPtr;

   *length = 0;

   loc1 = flagLoc( argc, argv, flag );
   if (badFlag) return( (int *)NULL );

   loc2 = nextFlag( argc, argv, loc1 );
   if (badFlag) return( (int *)NULL );

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       *length += cntToken( argv[loc] );
      }

   if (*length < 1)
      {
       *length = 0;
       return( (int *) NULL );
      }

   if (!(iArray = intAlloc( "listArg", *length )))
      {
       badFlag = 0;
       *length = -1;
       return( (int *)NULL );
      }

   iPtr = iArray;

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       if (loc < MAXARGS) argMap[loc] = '1';

       n = cntToken( argv[loc] );

       for( i = 1; i <= n; i++ )
          {
           *iPtr++ = getTokenI( i, argv[loc] );
          }
      }

   return( iArray );
}

/***/
/* rListArg: allocates and returns a float argument flag list;
/*           accomodates multiple items within each argument in argv.
/***/

float *rListArg( argc, argv, flag, length )

   int  argc, *length;
   char **argv, *flag;
{
   int   i, n, loc, loc1, loc2;
   float *rArray, *rPtr;

   *length = 0;

   loc1 = flagLoc( argc, argv, flag );
   if (badFlag) return( (float *)NULL );

   loc2 = nextFlag( argc, argv, loc1 );
   if (badFlag) return( (float *)NULL );

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       *length += cntToken( argv[loc] );
      }

   if (*length < 1)
      {
       *length = 0;
       return( (float *)NULL );
      }

   if (!(rArray = fltAlloc( "listArg", *length )))
      {
       badFlag = 0;
       *length = -1;
       return( (float *)NULL );
      }

   rPtr = rArray;

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       if (loc < MAXARGS) argMap[loc] = '1';

       n = cntToken( argv[loc] );

       for( i = 1; i <= n; i++ )
          {
           *rPtr++ = getTokenR( i, argv[loc] );
          }
      }

   return( rArray );
}

/***/
/* sListArg: allocates and returns a string argument flag list;
/*           accomodates multiple items within each argument in argv;
/*           free result with freeStrList().
/***/

char **sListArg( argc, argv, flag, length )

   int  argc, *length;
   char **argv, *flag;
{
   int   i, j, n, loc, loc1, loc2;
   char  **sArray, *sPtr;

   *length = 0;

   loc1 = flagLoc( argc, argv, flag );
   if (badFlag) return( (char **)NULL );

   loc2 = nextFlag( argc, argv, loc1 );
   if (badFlag) return( (char **)NULL );

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       *length += cntToken( argv[loc] );
      }

   if (*length < 1)
      {
       *length = 0;
       return( (char **)NULL );
      }

   sArray = (char **) voidAlloc( "strDup", sizeof(char *)*(*length + 1) );

   if (!sArray)
      {
       badFlag = 0;
       *length = -1;
       return( (char **)NULL );
      }

   j = 0;

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       if (loc < MAXARGS) argMap[loc] = '1';

       n = cntToken( argv[loc] );

       for( i = 1; i <= n; i++ )
          {
           sArray[j++] = strDup( getToken( i, argv[loc] ));
          }
      }

   sArray[*length] = (char *)NULL;

   return( sArray );
}

/***/
/* rListArg2: allocates and returns a float argument flag list from a single
/*            item in the command line.
/***/

float *rListArg2( argc, argv, flag, length )

   int  argc, *length;
   char **argv, *flag;
{
   int   i, loc;
   float *rArray, *rPtr;

   *length = 0;

   loc = flagLoc( argc, argv, flag );
   if (badFlag) return( (float *) NULL );
   if (loc > argc - 2) return( (float *) NULL );

   loc++;

   *length = cntToken( argv[loc] );

   if (loc < MAXARGS) argMap[loc] = '1';

   if (!(rArray = fltAlloc( "listArg", *length )))
      {
       badFlag = 0;
       *length = -1;
       return( (float *) NULL );
      }

   rPtr = rArray;

   for( i = 0; i < *length; i++ )
      {
       *rPtr++ = getTokenR( i+1, argv[loc] );
      }

   return( rArray );
}

/***/
/* charListArg: extract a string argument list; returns number of args in list.
/***/

int charListArg( argc, argv, flag, str, nMax )

   int  argc, nMax;
   char **argv, *flag, *str;
{
   int  n, nStr, nArgs, loc, loc1, loc2;
   char **sArray, **sPtr;

   loc1 = flagLoc( argc, argv, flag );
   if (badFlag) return( -1 );

   loc2 = nextFlag( argc, argv, loc1 );
   if (badFlag) return( -1 );

   *str  = '\0';
   nStr  = 0;
   nArgs = 0;

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       n = istrlen( argv[loc] );
       if (nArgs) n++;

       if (n + nStr >= nMax) return( nArgs );
 
       if (loc < MAXARGS) argMap[loc] = '1';

       if (nArgs) 
          {
           (void) strcat( str, " " );
           nStr++;
          }

       (void) strcat( str, argv[loc] );

       nStr += n;
       nArgs++;
      } 
       
   return( nArgs );
}

/***/
/* strListArg: allocates and returns a string argument flag list.
/***/

char **strListArg( argc, argv, flag, length )

   int  argc, *length;
   char **argv, *flag;
{
   int  loc, loc1, loc2;
   char **sArray, **sPtr;

   *length = 0;

   loc1 = flagLoc( argc, argv, flag );
   if (badFlag) return( (char **) NULL );

   loc2 = nextFlag( argc, argv, loc1 );
   if (badFlag) return( (char **) NULL );

   *length = loc2 - loc1 - 1;

   if (*length < 1)
      {
       *length = 0;
       return( (char **) NULL );
      }

   sArray = (char **) voidAlloc( "strListArg", sizeof(char *)*(*length + 1) );

   if (!sArray)
      {
       badFlag = 0;
       *length = -1;
       return( (char **) NULL );
      }

   sArray[*length] = (char *) NULL;
   sPtr = sArray;

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       if (loc < MAXARGS) argMap[loc] = '1';

       if (!(*sPtr = strDup( argv[loc] )))
          {
           FPR( stderr, "StrList Error allocating string.\n" );
          }

       sPtr++;
      } 
       
   return( sArray );
}

/***/
/* ptrListArg: returns original pointer to a string argument flag list.
/***/

char **ptrListArg( argc, argv, flag, length )

   int  argc, *length;
   char **argv, *flag;
{
   int  loc, loc1, loc2;
   char **sArray, **sPtr;

   *length = 0;

   loc1 = flagLoc( argc, argv, flag );
   if (badFlag) return( (char **) NULL );

   loc2 = nextFlag( argc, argv, loc1 );
   if (badFlag) return( (char **) NULL );

   *length = loc2 - loc1 - 1;

   for( loc = loc1 + 1; loc < loc2; loc++ )
      {
       if (loc < MAXARGS) argMap[loc] = '1';
      }

   return( argv + loc1 + 1 );
}

/***/
/* intArgD: sets an integer value if a flag is present;
/*          returns 1 if value was set.
/***/

int intArgD( argc, argv, flag, value )

   int  argc;
   char **argv, *flag; 
   int *value;
{
    int temp;

    temp = intArg( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        *value = temp;
        return( 1 );
       }
}

/***/
/* nmrIntArgD: sets an NMR_INT value if a flag is present;
/*             returns 1 if value was set.
/***/

int nmrIntArgD( argc, argv, flag, value )

   int  argc;
   char **argv, *flag;
   NMR_INT *value;
{
#ifdef NMR64
    char *temp;

    temp = ptrArg( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        *value = atol( temp );
        return( 1 );
       }
#else
    int temp;

    temp = intArg( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        *value = temp;
        return( 1 );
       }
#endif
}

/***/
/* logArgD: sets an integer logical value if a flag is present;
/*          returns 1 if value was set.
/***/

int logArgD( argc, argv, flag, value )

   int  argc;
   char **argv, *flag; 
   int *value;
{
    (void) flagLoc( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        *value = 1;
        return( 1 );
       }
}

/***/
/* notLogArgD: clears an integer logical value if a flag is present;
/*             returns 1 if value was cleared.
/***/

int notLogArgD( argc, argv, flag, value )

   int  argc;
   char **argv, *flag;
   int *value;
{
    (void) flagLoc( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        *value = 0;
        return( 1 );
       }
}

/***/
/* fltArgD: sets a float value if a flag is present;
/*          returns 1 if value was set.
/***/

int fltArgD( argc, argv, flag, value )

   int   argc;
   char  **argv, *flag; 
   float *value;
{
    float temp;

    temp = fltArg( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        *value = temp;
        return( 1 );
       }
}

/***/
/* dblArgD: sets a double value if a flag is present;
/*          returns 1 if value was set.
/***/

int dblArgD( argc, argv, flag, value )

   int    argc;
   char   **argv, *flag; 
   double *value;
{
    double temp;

    temp = dblArg( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        *value = temp;
        return( 1 );
       }
}

/***/
/* strArgD: copies a string value if a flag is present;
/*          returns 1 if value was set.
/***/

int strArgD( argc, argv, flag, value )

   int   argc;
   char  **argv, *flag; 
   char  *value;
{
    char *temp;

    temp = strArg( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        (void) strcpy( value, temp );
        return( 1 );
       }
}

/***/
/* ptrArgD: returns pointer to string value if a flag is present;
/*          returns 1 if value was set.
/***/

int ptrArgD( argc, argv, flag, value )

   int   argc;
   char  **argv, *flag;
   char  **value;
{
    char *temp;

    temp = ptrArg( argc, argv, flag );

    if (badFlag)
       {
        return( 0 );
       }
    else
       {
        *value = temp;
        return( 1 );
       }
}

/***/
/* parseLocList: set locations in iList to iON if they are listed in 
/*               command line as values or ranges from 1 to n, for 
/*               example: 2 5 7-10 15 ...
/***/

int parseLocList( argc, argv, flag, iList, n, iON, iOFF )

   int  argc, *iList, iON, iOFF, n;
   char *flag, **argv;
{
    int  i, ii, i1, iN, nON, nT, nN, iTemp, loc, loc1, locN;
    char *sPtr, *tPtr, *cPtr;

    if (!iList || n < 1) return( 0 );
   
    for( i = 0; i < n; i++ ) iList[i] = iOFF;

    loc1 = flagLoc( argc, argv, flag );
    if (loc1 < 0 || loc1 > argc - 2) return( 0 );

    locN = nextFlag( argc, argv, loc1 );
    nON  = 0;
 
    for( loc = loc1 + 1; loc < locN; loc++ )
       {
        sPtr = getNthArg( argc, argv, loc );
        nT   = cntToken( sPtr );

        for( i = 1; i <= nT; i++ )
           {
            tPtr = getToken( i, sPtr );
            nN   = cntTokenC( tPtr, "-" );

            if (nN == 1)
               {
                ii = 0;

                (void) sscanf( tPtr, "%d", &ii );
                if (ii >= 1 && ii <= n) iList[ii-1] = iON;
               }
            else if (nN > 1)
               {
                i1 = 0; iN = 0;

                cPtr = getTokenC( 1, tPtr, "-" );
                (void) sscanf( cPtr, "%d", &i1 );
                
                cPtr = getTokenC( 2, tPtr, "-" );
                (void) sscanf( cPtr, "%d", &iN );
         
                if (i1 > iN) {iTemp = i1; i1 = iN; iN = iTemp;} 

                for( ii = i1; ii <= iN; ii++ )
                   {
                    if (ii >= 1 && ii <= n) iList[ii-1] = iON;
                   }
               }
           }
       }
 
    for( i = 0; i < n; i++ ) 
       {
        if (iList[i] == iON)
           {
            nON++;
           }
       }

    return( nON );
}

/***/
/* copyArgv: returns a copy of argv, but with a new argv[0].
/***/

char **copyArgv( argc, argv, argv0 )

   char **argv, *argv0;
   int  argc;
{
    char **newArgv;
    int  i;

    if (!(newArgv = (char **)voidAlloc( "copyarg", sizeof(char *)*(argc+1) )))
       {
        return( (char **)NULL );
       }

    newArgv[0]    = strDup( argv0 );
    newArgv[argc] = (char *)NULL;

    for( i = 1; i < argc; i++ ) newArgv[i] = strDup( argv[i] ); 

    return( newArgv );
}

/***/
/* getFmtCount: returns number of single % characters in str.
/***/

int getFmtCount( str )

   char *str;
{
    int i, n, count;

    if (!str) return( 0 );

    n     = (int)strlen( str );
    count = 0;

    for( i = 0; i < n - 1; i++ )
       {
        if (str[i] != '%') continue;
        
        if (str[i+1] == '%')
           i++;
        else
           count++;
       }

   return( count );
}

/***/
/* isZFmt: returns 1 if template has only one % format specifier.
/***/

int isZFmt( template )

   char *template;
{
    int fmtCount;

    fmtCount = getFmtCount( template );

    if (fmtCount == 1) return( 1 );

    return( 0 );
}

/***/
/* isAZFmt: returns 1 if template has both A and Z format specifiers.
/***/

int isAZFmt( template )

   char *template;
{
    int fmtCount;

    fmtCount = getFmtCount( template );

    if (fmtCount == 2) return( 1 );

    return( 0 );
}


/***/
/* isBlank: return 1 if line is blank.
/***/

static int isBlank( line )

   char *line;
{
    if (!line) return( 1 );

    while( *line )
       {
        if (!isspace( *line )) return( 0 );
        line++;
       }

    return( 1 );
}

/***/
/* safe_strcpy: string copy, checking for null pointers.
/***/

char *safe_strcpy( dest, src )

   char *dest, *src;
{
    if (!dest) return( (char *)NULL );
  
    if (!src)
       dest[0] = '\0';
    else
       (void) strcpy( dest, src );

    return( dest );
}

/***/
/* envIsTrue: return 1 if environment variable is defined as 1, yes, or true.
/***/

int envIsTrue( eStr )

   char *eStr;
{
    char *vStr;

    if (!(vStr = getenv( eStr ))) return( 0 );

    if (!strcasecmp( vStr, "1" ))    return( 1 );
    if (!strcasecmp( vStr, "true" )) return( 1 );
    if (!strcasecmp( vStr, "yes" ))  return( 1 );

    return( 0 );
}

/***/
/* Bottom.
/***/
