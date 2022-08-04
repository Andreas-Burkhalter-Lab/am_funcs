/* son.h
**                                                               78 cols --->*
** Definitions of the structures and routines for the SON filing             *
** system. This is the include file for standard use, all access             *
** is by means of functions - no access to the internal SON data. The file   *
** machine.h provides some common definitions across separate platforms,     *
** note that the MS variants are designed to accomodate medium model - far   *
** pointers are explicitly used where necessary.
**
** SONAPI Don't declare this to give a pascal type on the Mac, there is a MPW
**        compiler bug that corrupts floats passed to pascal functions!!!!!!
**
*/

#ifndef __SON__
#define __SON__

#include "machine.h"

#if defined(macintosh) || defined(_MAC) /* define SONCONVERT in here if you want it */
    #include <Types.h>
    #include <Files.h>
    #include <Errors.h>
    #define USEHANDLES
#define  SONAPI(type) type
#undef LLIO                         /* LLIO is not used for Mac             */
#endif                              /* End of the Mac stuff, now for DOS    */

#ifdef _IS_MSDOS_
    #define qDebug 0                /* only used to debug Mac stuff         */
    #undef  USEHANDLES
    #include <malloc.h>
    #include <dos.h>
    #include <errno.h>
    #define LLIO                    /* We can use LLIO for MSC/DOS          */
    #define SONAPI(type) type _pascal
#endif

#if defined(_IS_WINDOWS_) && !defined(_MAC)
    #define qDebug 0                /* only used to debug Mac stuff         */
    #define USEHANDLES              /* always use handles under windows     */
    #define LLIO                    /* We can use LLIO for MSC/Windows      */
    #define SONAPI(type) type WINAPI
#endif

#define SONMAXCHANS 32      /* The app needs to know the limit on channels */

/*
** Now define the error constants used throughout the program, we start
** with the DOS errors returned by the MSDOS.ASM library
*/
#define SON_NO_FILE  -1
#define SON_NO_DOS_FILE -2
#define SON_NO_PATH -3
#define SON_NO_HANDLES -4
#define SON_NO_ACCESS  -5
#define SON_BAD_HANDLE -6
#define SON_MEMORY_ZAP -7
#define SON_OUT_OF_MEMORY -8
#define SON_INVALID_DRIVE -15   /* I think this is a Mac error TDB */
#define SON_OUT_OF_HANDLES -16  /* This refers to SON file handles ... */

/*
**  Nick added a new error code, used by OpenOldFile if the active file
**  is already open, this should not clash with any Mac OSerrs
*/
#define SON_FILE_ALREADY_OPEN -600

/*
** These two codes are not supplied by MSDOS.ASM, but stem from it's use,
** sadly the ReadHandle and WriteHandle routines return 0 for an error
*/
#define SON_BAD_READ -17
#define SON_BAD_WRITE -18

/*
** now some of our own errors, put in holes that we think we will never
** get from DOS... famous last words!
*/
#define SON_NO_CHANNEL -9
#define SON_CHANNEL_USED -10
#define SON_CHANNEL_UNUSED -11
#define SON_PAST_EOF -12
#define SON_WRONG_FILE -13
#define SON_NO_EXTRA -14
#define SON_CORRUPT_FILE -19
#define SON_PAST_SOF -20
#define SON_READ_ONLY -21

/*
** These constants define the number and length of various strings
*/
#define SON_NUMFILECOMMENTS 5
#define SON_COMMENTSZ 79
#define SON_CHANCOMSZ 71
#define SON_UNITSZ 5
#define SON_TITLESZ 9

/*
** These define the types of data we can store in our file, and a type
** that will hold one of these values.
*/
#define ChanOff 0
#define Adc 1
#define EventFall 2
#define EventRise 3
#define EventBoth 4
#define Marker 5
#define AdcMark 6       /* this is marker plus Adc data attached */
#define RealMark 7      /* Marker with real numbers attached */
#define TextMark 8      /* Marker with text attached */

typedef unsigned char TDataKind;      /* This holds one of the above values */

/*
** These constants defines the state of a created file
*/
#define FastWrite 0
#define NormalWrite 1

/*
**  The TMarker structure defines the marker data structure, which holds
**  a time value with associated data. The TAdcMark structure is a marker
**  with attached array of ADC data. TRealMark and TTextMark are very
**  similar - with real or text data attached.
*/
typedef long TSTime;
typedef short TAdc;
typedef char TMarkBytes[4];

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
#pragma pack(2)
#endif

typedef struct
{
    TSTime mark;        /* the marker times, as for events */
    TMarkBytes mvals;   /* the marker values */
} TMarker;


typedef struct
{
    TMarker m;          /* the marker structure */
    TAdc a[8000];       /* the attached ADC data */
} TAdcMark;


typedef struct
{
    TMarker m;          /* the marker structure */
    float r[80];        /* the attached floating point data */
} TRealMark;


typedef struct
{
    TMarker m;          /* the marker structure */
    char t[250];        /* the attached text data */
} TTextMark;

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
#pragma pack()
#endif

typedef TAdc FAR * TpAdc;
typedef TSTime FAR *TpSTime;
typedef TMarker FAR * TpMarker;
typedef TAdcMark FAR * TpAdcMark;
typedef char FAR * TpStr;
typedef const char FAR * TpCStr;
typedef WORD FAR * TpWORD;
typedef BOOLEAN FAR * TpBOOL;
typedef float FAR * TpFloat;


/* 
**
** The marker filter extensions to SON
**
** The declaration of the types is only to allow you to declare a structure
** of type TFilterMask. This structure is to be passed into the two
** functions only. No other use should be made, except to save it, if
** needed.
**
*/

typedef unsigned char TFilterElt;   /* element of a map */
typedef TFilterElt TLayerMask[32];  /* 256 bits in the bitmap */

typedef struct
{
    char cAllSet[4];        /* set non-zero if all markers enabled */
    TLayerMask aMask[4];    /* set of masks for each layer */
} TFilterMask;

typedef TFilterMask FAR *TpFilterMask;

#define SON_FALLLAYERS  -1

#define SON_FALLITEMS   -1

#define SON_FCLEAR      0
#define SON_FSET        1
#define SON_FINVERT     2
#define SON_FREAD       -1


#ifdef __cplusplus
extern "C" {
#endif


/*
** Now definitions of the functions defined in the code
*/
SONAPI(void) SONInitFiles( void );

#if defined(macintosh) || defined(_MAC)

SONAPI(short) SONOpenOldFile(ConstStr255Param name, short vRefNum, long dirID,
                    SignedByte perm);
SONAPI(short) SONOpenNewFile(ConstStr255Param name, short fMode, WORD extra,
                short vRefNum, long dirID, SignedByte perm, 
                OSType creator, OSType fileType);
#else

SONAPI(short) SONOpenOldFile (TpStr name, int iOpenMode);
SONAPI(short) SONOpenNewFile (TpStr name, short fMode, WORD extra) ;

#endif

SONAPI(BOOLEAN) SONCanWrite(short fh);
SONAPI(short) SONCloseFile (short fh);
SONAPI(short) SONEmptyFile(short fh);
SONAPI(short) SONSetBuffSpace(short fh);
SONAPI(short) SONGetFreeChan(short fh);
SONAPI(void) SONSetFileClock(short fh, WORD usPerTime, WORD timePerADC);
SONAPI(short) SONSetADCChan(short fh, WORD chan, short adcChan, short dvd,
                 short buffSz, TpCStr com, TpCStr title, float ideal,
                 float scl, float offs, TpCStr unt);
SONAPI(short) SONSetADCMarkChan(short fh, WORD chan, short adcChan, short dvd,
                 short buffSz, TpCStr com, TpCStr title, float ideal, float scl,
                 float offs, TpCStr unt, WORD points, short preTrig);
SONAPI(short) SONSetRealMarkChan(short fh, WORD chan, short phyChan,
                 short buffSz, TpCStr com, TpCStr title, float ideal,
                 float min, float max, TpCStr unt, WORD points);
SONAPI(short) SONSetTextMarkChan(short fh, WORD chan, short phyChan,
                 short buffSz, TpCStr com, TpCStr title,
                 float ideal, TpCStr unt, WORD points);
SONAPI(void) SONSetInitLow(short fh, WORD chan, BOOLEAN bLow);
SONAPI(short) SONSetEventChan(short fh, WORD chan, short evtChan,
                 short buffSz, TpCStr com, TpCStr title, float ideal, TDataKind evtKind);
SONAPI(short) SONUpdateStart(short fh);
SONAPI(void) SONSetFileComment(short fh,WORD which, TpCStr comment);
SONAPI(void) SONGetFileComment(short fh,WORD which, TpStr comment, short max);
SONAPI(void) SONSetChanComment(short fh, WORD chan, TpCStr comment);
SONAPI(void) SONGetChanComment(short fh,WORD chan, TpStr comment, short max);
SONAPI(void) SONSetChanTitle(short fh, WORD chan, TpCStr title);
SONAPI(void) SONGetChanTitle(short fh,WORD chan, TpStr title);
SONAPI(void) SONGetIdealLimits(short fh, WORD chan, TpFloat ideal, TpFloat min,
                 TpFloat max);
SONAPI(WORD) SONGetusPerTime(short fh);
SONAPI(WORD) SONGetTimePerADC(short fh);
SONAPI(void) SONSetADCUnits(short fh, WORD chan, TpCStr units);
SONAPI(void) SONSetADCOffset(short fh, WORD chan, float offset);
SONAPI(void) SONSetADCScale(short fh, WORD chan, float scale);
SONAPI(void) SONGetADCInfo(short fh, WORD chan, TpFloat scale, TpFloat offset,
                 TpStr units, TpWORD points, short FAR *preTrig);
SONAPI(void) SONGetExtMarkInfo(short fh, WORD chan, TpStr units,
                 TpWORD points, short FAR *preTrig);
SONAPI(short) SONWriteEventBlock(short fh,WORD chan, TpSTime buffer,
                        long count);
SONAPI(short) SONWriteMarkBlock(short fh, WORD chan, TpMarker buffer,
                        long count);
SONAPI(TSTime) SONWriteADCBlock(short fh, WORD chan, TpAdc buffer,
                        long count, TSTime sTime);
SONAPI(short) SONWriteExtMarkBlock(short fh, WORD chan, TpMarker buffer,
                        long count);
SONAPI(short) SONCommitFile(short fh);
SONAPI(long) SONGetEventData(short fh,WORD chan,TpSTime evtDataP, long max,
                  TSTime sTime, TSTime eTime, TpBOOL levLowP, 
                  TpFilterMask pFiltMask);
SONAPI(long) SONGetMarkData(short fh,WORD chan,TpMarker markP, long max,
                  TSTime sTime,TSTime eTime, TpFilterMask pFiltMask);
SONAPI(long) SONGetADCData(short fh,WORD chan,TpAdc adcDataP, long max,
                  TSTime sTime,TSTime eTime,TpSTime pbTime, 
                  TpFilterMask pFiltMask);
SONAPI(long) SONGetExtMarkData(short fh,WORD chan, TpMarker markP, long max,
                  TSTime sTime,TSTime eTime, TpFilterMask pFiltMask);
SONAPI(long) SONGetExtraDataSize(short fh);
SONAPI(int) SONGetVersion(short fh);
SONAPI(short) SONGetExtraData(short fh, void FAR *buff, WORD bytes,
                  WORD offset, BOOLEAN writeIt);
SONAPI(short) SONSetMarker(short fh, WORD chan, TSTime time, TpMarker markP,
                  WORD size);
SONAPI(short)  SONChanDelete(short fh, WORD chan);
SONAPI(TDataKind) SONChanKind(short fh, WORD chan);
SONAPI(TSTime) SONChanDivide(short fh, WORD chan);
SONAPI(WORD)   SONItemSize(short fh, WORD chan);
SONAPI(short)  SONGetNewFileNum (void);
SONAPI(TSTime) SONChanMaxTime(short fh, WORD chan);
SONAPI(TSTime) SONMaxTime(short fh);

SONAPI(TSTime) SONLastTime(short fh, WORD channel, TSTime sTime, TSTime eTime,
                    short FAR * psVal, TMarkBytes FAR *mBytes,
                    TpBOOL pbMk, TpFilterMask pFiltMask);


SONAPI(int) SONFilter(TpMarker pM, TpFilterMask pFM);
SONAPI(int) SONFControl(TpFilterMask pFM, int layer, int item, int set);
SONAPI(BOOLEAN) SONFEqual(TpFilterMask pFiltMask1, TpFilterMask pFiltMask2);

#ifdef __cplusplus
}
#endif

#endif /* __SON__ */
