/*
** son.c
**
******************************************************************************
**
**  Son filing section main code
**
**  Copyright (c) Cambridge Electronic Design Limited 1988,1990
**
**  Author Greg P Smith
**
**  Revision history:
**
** 13/Sep/90 GPS The first C version of the library seems to work but not
**               yet tested. It emulates the Pascal version and further
**               isolates the user from the internal data strucures.
**
** 30/Apr/91 NJC the file is now being converted on the Mac.
**               I could also include code for doing lower level calls with
**               PBRead etc. via param blocks, and async, sync routines
**               The Mac version returns the value of the Mac OSerr code, not
**               the MSDOS error.
**               NB Always use short and long (MAC int is 32 bits, IBM is 16)
**               There are a few filey things on the Mac which may affect
**               future versions of the library - shared filing systems, and
**               sys 7 stuff for live sharing of data between different
**               applications, Eg if Spike2 data capture program is running
**               concurrently with analysis program. This will need more work
**               interface may change.
**
**          DIFFERENCES between DOS and Mac
**
**               The main differences between DOS and MAC versions are controlled
**               by the machine.h file. This is included in SON.H.
**
**               file handle field has been renamed refNum,(there is a type
**               called Handle on Mac) NB. In all routines which set file
**               refNum, you must pass short *, as they are declared in pascal
**               as VAR integer, ie 16 bits
**
**               FAR, _near are all defined as nothing if macintosh
**
**               file names are still passed as char *, however they must be
**               pascal strings with a length byte to start with. The calling
**               routines must ensure this is the case
**
**          POSSIBLE PROBLEM on DOS/Windows
**               SONWriteBlock() makes use of a global called outChannelP, by
**               means of an #ifdef SONCONVERT. On Windows we may have to
**               change that. This relies on the variable being correctly
**               set up by the caller.
**
**          VERSION NO. OF FILING SYSTEM and when to swap.
**               This needs properly specifying before first verison of
**               library ships. Currently, version 3 DOS files have garbage
**               in the osFormat word, which I have added. (it used to be
**               padding.) If we change to version 4 we should ensure that
**               the word is either 0x0000 (DOS) or 0x0101 (MAC)
**
**               When swapping a float value, it is first necessary to ensure
**               it is 4 byte aligned. This is done here by the procedure
**               FloatSwap which takes a float passed by value (hence copied
**               to a 4 byte boundary)
**
** 06/Jun/91 NJC Original SONWriteBlock contained a bug with & instead of &&
**
** 02/Aug/91     SONFindDataChan 2 bugs fixed
**
** 04/Jul/91     Discovered bug in SONGetSucc, SONGetPred. They can
**               return SON_NO_FILE (ie -1) or -1 for end of chain. Requires that
**               SON_NO_FILE be redefined. We could make a return of 0 from
**               SONGetSucc etc be an error
**
** 25/10/91      SONGetBuffSpace bug- checked phySz of all chans except 0.
**
**               SONGetNewFileNum () added, which returns the next available
**               file index number, (ie a no. in range 0..MAXFILES-1). The
**               application does not then need to keep track of which indices
**               have been allocated - it simply calls this procedure.
**               Used especially by MacAp
**
** 8.8.91 Mark   added SONGetIdealRate
**
**               SONOpenOldFile does not call SONCloseActiveFile at the start.
**               Instead it checks if the active file is open, if so it returns
**               SON_FILE_ALREADY_OPEN error. SONOpenNewFile will fail with the
**               same error
**
**               Longswap for sonconvert is no longer static as used in
**               sonconvert.c
**
** 04/Mar/92 TDB Code now back on the PC - first of all tidied up for the PC
**               by fixing up all of the memory pointers so that the code
**               compiles under the small model, but uses far pointers
**               extensively. The various defines have been grouped together
**               and documented - for ease of use.
**
** 13/Mar/92 TDB Now uses file handle to reference file, activeFile
**               disappears at last! Brought up to version 4 level, which
**               provides ADC marker data. FiltMark is now a procedural
**               parameter, not a global variable. REVISION = 4.
**
** 19/Mar/92 TDB Now extended still further to provide Real marker and
**               text marker data types with appropriate extra information.
**               Read and write routines for all extended marker data types
**               are now in common. REVISION = 5. More far parameters.
**
** 26/Jun/92 GPS Tidied up somewhat. Fixed silly in SONCommit...
**
** 10/Aug/92 GPS Changes for SONCONVERT from Mark.
**
**  8/Oct/92 MAE Changed SpeedPtrs to be allocated dynamically for each file,
**               and checked this change out for Windows.
**
**  3/Nov/92 MAE Added functions SONMaxTime and SONChanMaxTime (moved from
**               SonUtils.c)
**
**  4/Nov/92 PC  Changed so that max times are actually stored in the son file
**               (channel comment size decreased by 8 to provide space).
**
**  5/Nov/92 NJC Global var convert (only used by SonConvert) renamed
**               gConvert.Also tidied up the printfs in SonConvert.
**
**  4/Dec/92 MAE Added file handle to parameter list of FMarker filter 
**               function type definition.
**
** 24/Feb/93 PNC Changed to use new smybols defined in MACHINE.H for
**               different machine/version builds - WINNT now supported.
**
**  8/Mar/93 NJC Added global var gOrigVersion, (used only by SonConvert)
**               which is set to the file version, before upgrading takes
**               place
**
**  6/Apr/93 MAE SONSetBuffSpace only reallocates mem if the buffer size
**               has changed.
**               Also added function SONEmptyFile to remove all data from an 
**               existing son file.
**
** 28/Apr/93 MAE Fixed bug in SONFindBlock, which could return a block 
**               (instead of 0) if you asked for data entirely before the
**               first block.
**
** 22/Jun/93 MAE Added SONExtendMaxTime and SONGetFirstData for Windows,
**               plus SONFileRefNum so global 'files' does not need to
**               be accessed.
**
**  6/Jul/93 MAE SONGetEventData was returning the wrong level (the opposite 
**               of what it should) except when asking for data before the
**               first block.
**
** 23/Jul/93 TDB Adjusted many functions to ensure parameters not used before
**               ASSERTs. Removed redundant error codes. NJC fixed Mac form
**               of SONCommitFile. GetNewFileNum, SONReOpenOldFile are now
**               just internal routines. SONFindDataChan removed. Use of
**               long instead of TSTime corrected. SONGetFileClock removed,
**               SONGetTimePerADC added. SONSetChanComment, SONSetInitLow
**               added. float* changed to float FAR * (TpFloat), ditto
**               FMARKER* changed to defined far pointer type TpFMarker.
**               Experimental APPENDWRITE code added (not tested) and new
**               function SONSetMarker written.
**
** 28/Sep/93 MAE Added in Gregs code for filtering of markers. Functions now 
**               take a pointer to a filter mask structure (TpFilterMask) 
**               instead of a filter function.
**
**  4/Oct/93 MAE Fixed bug in SONGetAdcData when called for an AdcMark 
**               channel. If the end time was within the adc data of the
**               marker then you were returned one point too few.
**
** 26/Oct/93 PNC Fixed bug causing unnecessary writes; updateHead was
**               called in Close for non-SON files and for files internally
**               updated to V5 with no other writes needed.
**
** 14/Dec/93 MAE Added SONSetChanTitle & SONSetADCUnits. Also made sure that
**               these 2 and the following set updateHead: SONSetFileComment,
**               SONSetChanComment, SONSetFileClock, SONSetInitLow.
**               Also modified SONFControl so that SON_FALLLAYERS will work
**               without SON_FALLITEMS as expected.
**
** 25/Jan/94 MAE Added end of search time to SONLastTime (note that this is
**               earlier than the start time).
**
**  4/May/94 MAE Fixed bug in SONChanDelete - if you deleted a channel with
**               no data blocks that did have some deleted blocks then the
**               deleted chain pointers would be corrupted. Converted a pair
**               of near pointers in filter manipulation to FAR to avoid
**               segment loss in WIN16 build.
**
**  2/Dec/94 TDB Fixed fault in all SONWriteXXXXBlock routines; variable that
**               indexed through data was WORD while count of items was long.
**
**  2/Dec/94 MAE Added SONSetADCOffset and SONSetADCScale.
**               Fix in SONFindBlock so won't return a block if end < start.
**
** 22/Dec/94 MAE Added check in SONSetMarker for succBlock <= lastBlock for
**               the case when this is called during sampling and the next
**               block has not been writen yet.
**
** 22/Mar/95 NJC Added in an error return from ConvertDataBlock,
**               called in SONReadBlock (SonConvert only) 
**               Removed the call to SONUpdateMaxTimes (SonConvert only)
**               called in SONReOpenOldFile to upgrade from sys 4 to sys 5
**               when converting from DOS to Mac, and the call added to
**               SONCloseFile
**               SONCloseFile now calls FlushVol whether or not there was
**               a prior disk error (macintosh only)
**
**  1/Sep/95 TDB Altered SONINTL.H to increase max files to 32.
**
**  8/Nov/95 TDB Fixed SonFilter to avoid problems with marker values
**               that were greater than 127 - were treated as negative.
**
** 14/Oct/96 TDB Altered variable called this, and added casts to TMarker*,
**               both to get SON to compile when compiler is in c++ mode.
**
** 14/Jul/98 GPS Added SONGetExtraDataSize(short fh) to get size of the extra
**               data area, and SONGetVersion(short fh) to return the file
**               version number. SONCanWrite returns TRUE if file is rd/write.
**
*/

#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <limits.h>

/*************************************************************************
**
** These are the machine specific definitions and includes. Some of them
** are not really machine specific, but I've put them here anyway. first
** of all are the system defines, next machine global, finally specific
**
** Remove the NDEBUG definition to cause the asserts to be used
** The LLIO (low level i/o) should be undefined to use streams library and
** defined for low level. This is to make it a bit easier to implement the
** code on different machines which may not have the low level stuff.
** LLIO is defined in son.h and used in son.h and here
**
**
** macintosh and MSC are set by the compiler, and define machine type
** _IS_MSDOS_     set in machine.h for an msdos native mode build
** _IS_WINDOWS_   set in machine.h for windows 16 or 32 bit mode
** qDebug         (mac)set by MPW to enable debug fprintf
**
** NDEBUG         gets rid of the asserts.
** USEHANDLES     To use memory handles in place of pointers. This form is
**                supported by code generated for the Mac and Windows forms.
**
** LLIO           (dos) if defined, low level I/O is used, otherwise streams.
** SONCONVERT     (mac) define if we should convert format (little-big endian)
**  [convert      (mac) specifies direction of conversion - for SonConvert]
**
*/

/*#undef   NDEBUG    No asserts if this is defined */
#if qDebug
#undef NDEBUG
#else
#define NDEBUG 1
#endif

#include <assert.h>
#include "Son.h"

#ifdef LLIO
#include <io.h>
#include <sys\types.h>
#include <sys\stat.h>
#include <fcntl.h>
#include <errno.h>
#endif

#include "SonIntl.h"

#ifdef SONCONVERT
    #include "Convert.h"
#endif

#ifdef macintosh
    #pragma segment Son
#endif

/******************* fixed data structures *****************************
** 
** These are the fixed data structures on which the library is based.
** workP points at an area used to read the header information from a
** data block for fast disk searching and filling in the forward and
** backward links.
*/

TSonFile _near files[MAXFILES];   /* the possible files we allow */
char _near workBlock[DISKBLOCK];  /* array of DISKBLOCK bytes for workP */
TpDataBlock _near workP = (TpDataBlock)workBlock; /* SONGetSucc, SONSetSucc */

#ifdef USEHANDLES
/**************************************************************************
**
** Handle support is for systems which use memory handles (Mac and Windows)
** and not simply pointers. The code assume Windows if USEHANDLES and
** not macintosh.
*/

/********************** H A l l o c ( ) **********************************
**
** Allocates space of the required size, sets up a handle and returns a
** pointer
*/
static void FAR *HAlloc(long size, THandle FAR *handleP)
{
    THandle hdl ;
#ifdef macintosh
    hdl = (THandle) NewHandle (size) ;
    if (hdl)
    {
        MoveHHi (hdl);
        HLock (hdl) ;
        *handleP = hdl ;
        return *hdl ;
#endif
#ifdef _IS_WINDOWS_
    hdl = GlobalAlloc (GMEM_MOVEABLE, size) ;
    if (hdl)
    {
        *handleP = hdl ;
        return GlobalLock (hdl) ;
#endif
    }
    else
        return NULL ;
}

/********************** H F r e e ( ) ************************************
**
** Frees space pointed to by theHandle
*/
static void HFree(THandle theHandle)
{
#ifdef macintosh
    DisposeHandle(theHandle);
#endif
#ifdef _IS_WINDOWS_
    GlobalUnlock(theHandle);
    GlobalFree(theHandle);
#endif
}
#endif  /* of USEHANDLES only bit */


/**************** S O N R e a d ( ) ************************************
**
** This is the generic son system function to read data from the current
** file to a buffer, from a given position, and to check for an error.
** The function returns 0 if all is ok, otherwise it returns an error.
**
** On the Mac, if buffer is a dereferenced handle, then handle must be
** locked. The Mac version returns a Mac OSErr code. Normally only used
** internally.
**
** fh       The handle for the file
** buffer   Pointer to the area data to be transferred to
** bytes    Number of bytes to be read/written
** offset   Position in the file to be used
*/
SONAPI(short) SONRead(short fh, void FAR * buffer, WORD bytes, long offset)
{
#if defined(macintosh) || defined(_MAC)
    OSErr err ;
    long nBytes = bytes ;               /* Mac read routine expects longint */
#endif

#ifdef _IS_MSDOS_
#ifdef LLIO
    unsigned rbytes ;
#endif
#endif

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(offset >= 0);

#if defined(macintosh) || defined(_MAC)
    err = SetFPos (files[fh].refNum, fsFromStart, offset) ;   /* seek */
    if (err)
    {
#if qDebug
        fprintf (stderr,"ERROR: %d SetFPos, offset %d\n",err,offset);
#endif
        return err ;
    }

    err = FSRead (files[fh].refNum, &nBytes, buffer) ;
    if (err)
    {
#if qDebug
        fprintf (stderr, "ERROR: %d FSRead, nbytes read %d\n", err, nBytes);
#endif
        return err ;
    }
#endif

#if defined(_IS_WINDOWS_) && !defined(_MAC)
    if (_llseek (files[fh].handle, offset, SEEK_SET) != offset)/* Windows read */
        return SON_BAD_READ;
    if (_lread(files[fh].handle, buffer, bytes) != bytes)
        return SON_BAD_READ;
#endif  /* if Windows */

#ifdef _IS_MSDOS_
#ifdef LLIO
    if (lseek (files[fh].handle, offset, SEEK_SET) != offset)/* LLIO read */
        return SON_BAD_READ;
    if (_dos_read(files[fh].handle, buffer, bytes, &rbytes) != 0)
        return SON_BAD_READ;
#else
    if (fseek (files[fh].handle, offset, SEEK_SET))     /* stdio read */
        return SON_BAD_READ;
    if (fread(buffer, 1, bytes, files[fh].handle) != bytes)
        return SON_BAD_READ;
#endif  /* if LLIO else */
#endif  /* if MSDOS else */
    return 0;
}

/****************** S O N W r i t e ( ) *********************************
**
** This is the generic son system function to write data from the current
** file to a buffer, from a given position, and to check for an error.
** The function returns 0 if all is ok, otherwise it returns an error.
**
** Normally only used internally.
**
** fh       The handle for the file
** buffer   Pointer to the area to be transferred
** bytes    Number of bytes to be read/written
** offset   Position in the file to be used
*/
SONAPI(short) SONWrite(short fh, void FAR * buffer, WORD bytes, long offset)
{
#if defined(macintosh) || defined(_MAC)
    OSErr err ;
    long nBytes = bytes;    /* Mac read routine expects longint */
#endif

#ifdef _IS_MSDOS_
#ifdef LLIO
    unsigned wbytes ;
#endif
#endif

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(offset >= 0);

    if (files[fh].bReadOnly)
        return SON_READ_ONLY;

#if defined(macintosh) || defined(_MAC)
    err = SetFPos (files[fh].refNum, fsFromStart, offset) ;
    if (err)
        return err ;
    err = FSWrite (files[fh].refNum, &nBytes, buffer) ;
    if (err)
        return err ;
#endif

#if defined(_IS_WINDOWS_) && !defined(_MAC)
    if (_llseek (files[fh].handle, offset, SEEK_SET) != offset)
        return SON_BAD_WRITE;
    if (_lwrite(files[fh].handle, (LPCSTR)buffer, bytes) != bytes)
        return SON_BAD_WRITE;
#endif /* if Windows */

#ifdef _IS_MSDOS_
#ifdef LLIO
    if (lseek (files[fh].handle, offset, SEEK_SET) != offset)
        return SON_BAD_WRITE;
    if (_dos_write(files[fh].handle, buffer, bytes, &wbytes) != 0)
        return SON_BAD_WRITE;
#else
    if (fseek (files[fh].handle, offset, SEEK_SET))
        return SON_BAD_WRITE;
    if (fwrite(buffer, 1, bytes, files[fh].handle) != bytes)
        return SON_BAD_WRITE;
#endif /* if LLIO else */
#endif /* if MSDOS */
    return 0;
}

/****************** s t r 2 l s t r ( ) ***********************************
**
** convert a string to a pascal type LSTRING, copying at most max chars.
** The receiving string (l) needs to be max+1 characters long
*/
static void str2lstr(void FAR *l, const char * s, int max)
{
    int len = F_strlen(s);      /* how long is the source */
    TpStr t = (TpStr) l;        /* local copy of pointer */
    
    if (len>max) len = max;     /* make sure not too long */
    *t++ = (char)len;           /* set the length */
    for (;len;len--)
        *t++ = *s++;            /* copy the string */
}

/***************** l s t r 2 s t r ( ) ************************************
**
** convert a pascal string to a string, copying at most max chars. max does
** not include the \0 on the end. Changed 8.8.91 to not include the \0 which
** makes it more compatible with str2lstr. The receiving string (s) needs to
** be max+1 characters long
*/
static void lstr2str(TpStr s, void FAR *l, int max)
{
    int len;    
    TpStr t = (TpStr) l;        /* local copy of pointer */
    
    len = *t++;                 /* get length, point at first character */
    if (len>max)
        len = max;              /* make sure not too long */
    for (;len;len--)
        *s++ = *t++;
    *s = '\0';                  /* terminate the string */
}


/***************** S O N I t e m S i z e ( ) **************************
**
** Returns size of channel data, in bytes, particularly useful with
** the extended marker types.
**
*/
SONAPI(WORD) SONItemSize(short fh, WORD chan)
{
    TpChannel pChan;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    pChan = SONChanPnt(fh,chan);                    /* point at our channel */
    switch (pChan->kind)
    {
    case Adc:
        return sizeof(TAdc);

    case Marker :
        return sizeof(TMarker);

    case AdcMark :
    case RealMark :
    case TextMark :
        return (WORD)(sizeof(TMarker) + pChan->nExtra);

    case EventFall :
    case EventRise :
    case EventBoth :
        return sizeof(TSTime);

    default : return 1;     /* to avoid divide by zero! */
    }
}

/***************** S O N C h a n P n t ( ) *******************************
**
** Used to get a pointer to an active channel. This saves quite a bit
** of code each time it is used.
*/
SONAPI(TpChannel) SONChanPnt(short fh, WORD chan)
{
    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    return files[fh].chanP+chan;           /* assume file and chan OK */
}

/***************** S O N C h a n K i n d ( ) ****************************
**
** A quick way to find out the type of a channel.
*/
SONAPI(TDataKind) SONChanKind(short fh, WORD chan)
{
    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    return SONChanPnt(fh, chan)->kind;
}

/***************** S O N S e t F i l e C o m m e n t ( ) ***************
**
** Write the comment associated with the file. Up to 5 lines of comment
** can be set (0-4 as the index).
*/
SONAPI(void) SONSetFileComment(short fh, WORD which, TpCStr comment)
{
    assert(which<SON_NUMFILECOMMENTS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (which < SON_NUMFILECOMMENTS)
    {
        str2lstr((void FAR *)&(files[fh].headP->fileComment[which]),
                                                 comment, SON_COMMENTSZ);
        files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
    }
}

/************** S O N G e t F i l e C o m m e n t ( ) ********************
**
** return the file comment in the string pointed at by comment. Null comment
** returned if <which> is silly.
*/
SONAPI(void) SONGetFileComment(short fh, WORD which, TpStr comment, short max)
{
    assert(which<SON_NUMFILECOMMENTS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (which < SON_NUMFILECOMMENTS)
        lstr2str(comment,(void FAR *)&(files[fh].headP->fileComment[which]),
                                                                         max);
    else
        *comment = (char)0;
}

/************** S O N S e t C h a n C o m m e n t ( ) ****************
**
** Set the comment for a channel. Provided separately so that comment can
** be changed at the end of sampling run, not fixed at the beginning.
*/
SONAPI(void) SONSetChanComment(short fh, WORD chan, TpCStr comment)
{
    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    str2lstr((void FAR *)&(SONChanPnt(fh, chan)->comment), comment, SON_CHANCOMSZ);
    files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
}

/************** S O N G e t C h a n C o m m e n t ( ) ****************
**
** Get the comment back from a channel
*/
SONAPI(void) SONGetChanComment(short fh, WORD chan, TpStr comment, short max)
{
    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    lstr2str(comment, (void FAR *)&(SONChanPnt(fh, chan)->comment), max);
}

/************** S O N S e t C h a n T i t l e ( ) ****************
**
** Set the title for a channel. Provided separately so that title can be
** changed once channel has been created.
*/
SONAPI(void) SONSetChanTitle(short fh, WORD chan, TpCStr title)
{
    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    str2lstr((void FAR *)&(SONChanPnt(fh, chan)->title), title, SON_TITLESZ);
    files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
}

/*************** S O N G e t C h a n T i t l e ( ) *******************
**
** Get the title for a channel
*/
SONAPI(void) SONGetChanTitle(short fh, WORD chan, TpStr title)
{
    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    lstr2str(title, (void FAR *)&(SONChanPnt(fh, chan)->title), SON_TITLESZ);
}

/*************** S O N G e t I d e a l L i m i t s ( ) ****************
**
** Get the ideal rate and expected Y axis limits, as appropriate;
*/
SONAPI(void) SONGetIdealLimits(short fh, WORD chan, TpFloat ideal,
                               TpFloat min, TpFloat max)
{
    TpChannel pC;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    pC = SONChanPnt(fh,chan);                    /* point at our channel */
    *ideal = pC->idealRate;                            /* Get ideal rate */
    switch (pC->kind)
    {
        case Adc:
        case AdcMark:               /* limits to ADC scaling */
        {
            *min = -5 * pC->v.adc.scale + pC->v.adc.offset;
            *max =  5 * pC->v.adc.scale + pC->v.adc.offset;
            break;
        }
        case EventFall:
        case EventRise:
        case EventBoth:
        case Marker:
        case TextMark:
        {
            *min = 0.0f;                   /* return expected limits to rate */
            *max = pC->idealRate;
            break;
        }
        case RealMark:
        {
            *min = pC->v.real.min;         /* return stored expected limits */
            *max = pC->v.real.max;
            break;
        }
        default :
        {
            *min = 0.0f;
            *max = 1.0f;
            *ideal = 1.0f;
            break;
        }
    }
}

/**************** S O N G e t u s P e r T i m e ( ) *******************
**
** Get usPerTime
*/
SONAPI(WORD) SONGetusPerTime (short fh)
{
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (!files[fh].defined)
        return 0;
    return files[fh].headP->usPerTime;
}

/*************** S O N G e t T i m e P e r A D C ( ) ******************
**
** Get TimePerADC
*/
SONAPI(WORD) SONGetTimePerADC (short fh)
{
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (!files[fh].defined)
        return 0;
    return files[fh].headP->timePerADC;
}

/************** S O N S e t A D C U n i t s ( ) ****************
**
** Set the units for an adc channel. Provided separately so that units
** string can be changed once channel has been created.
*/
SONAPI(void) SONSetADCUnits(short fh, WORD chan, TpCStr units)
{
    TpChannel cP;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    cP = SONChanPnt(fh, chan);                    /* Get pointer to channel */
    if ((cP->kind != Adc) && (cP->kind != AdcMark))
        return;                         /* Check that channel kind is valid */

    str2lstr((void FAR *)&(cP->v.adc.units), units, SON_UNITSZ);
    files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
}

/************** S O N S e t A D C O f f s e t ( ) ****************
**
** Set the offset for an adc channel. Provided separately so that offset
** can be changed once channel has been created.
*/
SONAPI(void) SONSetADCOffset(short fh, WORD chan, float offset)
{
    TpChannel cP;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    cP = SONChanPnt(fh, chan);                    /* Get pointer to channel */
    if ((cP->kind != Adc) && (cP->kind != AdcMark))
        return;                         /* Check that channel kind is valid */

    cP->v.adc.offset = offset;
    files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
}

/************** S O N S e t A D C S c a l e ( ) ****************
**
** Set the scale for an adc channel. Provided separately so that scale
** can be changed once channel has been created.
*/
SONAPI(void) SONSetADCScale(short fh, WORD chan, float scale)
{
    TpChannel cP;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    cP = SONChanPnt(fh, chan);                    /* Get pointer to channel */
    if ((cP->kind != Adc) && (cP->kind != AdcMark))
        return;                         /* Check that channel kind is valid */

    cP->v.adc.scale = scale;
    files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
}

/****************** S O N G e t A D C I n f o ( ) ********************
**
** Get the info for an ADC channel
*/
SONAPI(void) SONGetADCInfo(short fh, WORD chan, TpFloat scale, TpFloat offset,
                            TpStr units, TpWORD points, short FAR *preTrig)
{
    TpChannel cP;

    assert(chan < SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(scale != NULL);
    assert(offset != NULL);
    assert(units != NULL);
    assert(points != NULL);
    assert(preTrig != NULL);

    cP = SONChanPnt(fh, chan);                    /* Get pointer to channel */
    if ((cP->kind != Adc) && (cP->kind != AdcMark))
        return;                         /* Check that channel kind is valid */
    *scale = cP->v.adc.scale;
    *offset = cP->v.adc.offset;
    lstr2str(units, (void FAR *)&(cP->v.adc.units), SON_UNITSZ);
    *points = (WORD)((cP->kind == AdcMark) ? (cP->nExtra >> 1) : 1);
    *preTrig = (short)((cP->kind == AdcMark) ? cP->preTrig : 0);
}

/**************** S O N G e t E x t M a r k I n f o ( ) ***************
**
** Get the info for an extended marker channel. This will return information
** appropriate to all extended marker channels.
*/
SONAPI(void) SONGetExtMarkInfo(short fh, WORD chan, TpStr units,
                                       TpWORD points, short FAR *preTrig)
{
    TpChannel cP;

    assert(chan < SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(units != NULL);
    assert(points != NULL);
    assert(preTrig != NULL);

    cP = SONChanPnt(fh, chan);
    if ((cP->kind != RealMark) &&       /* Check that channel kind is valid */
       (cP->kind != AdcMark)   &&
       (cP->kind != TextMark))
        return;
    lstr2str(units, (void FAR *)&(cP->v.adc.units), SON_UNITSZ);
    switch (cP->kind)
    {
        case AdcMark :
        {
            *points = (WORD)(cP->nExtra >> 1);
            *preTrig =  cP->preTrig;
            break ;
        }
        case RealMark :
        {
            *points = (WORD)(cP->nExtra >> 2);
            *preTrig =  0;
            break ;
        }
        case TextMark :
        {
            *points = cP->nExtra;
            *preTrig =  0;
            break ;
        }
        default: ;
    }
}

/***************** S O N C h a n D i v i d e ( ) *********************
** SONChanDivide:- Get the divide down from the tick rate for an ADC channel
** This is returned as 1 if we ask for an event or off channel.
*/
SONAPI(TSTime) SONChanDivide(short fh, WORD chan)
{
    TDataKind nKind;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    nKind = SONChanKind(fh,chan);
    if ((nKind!=Adc) && (nKind!=AdcMark))   /* not an adc channel? */
        return 1;
    else
        return ((TSTime)files[fh].headP->timePerADC) *
                SONChanPnt(fh, chan)->v.adc.divide;
}

/***************** S O N I n t l C h a n M a x T i m e ( ) ******************
** SONIntlChanMaxTime:-  returns the end time contained in the last data
** block of the channel BY LOOKING THROUGH THE DATA ITSELF
** Return 0 if channel is off or deleted
*/
SONAPI(TSTime) SONIntlChanMaxTime(short fh, WORD chan)
{
    short err;
    long last;                  /* address of last block */
    TpChannel chanP;            /* pointer to channel block */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    chanP = SONChanPnt(fh, chan);
    if ( (SONChanKind (fh, chan) == ChanOff) )
        return 0;
    else                                /* get last block */
    {
        last = chanP->lastBlock;
        if (last == CHAINEND)           /* or end of chain */
            return 0;
        else
        {
            err = SONGetBlock (fh, last);
            if (err == 0)
                return workP->endTime;  /* we have found it */
            else
                return err;
        }
    }
}

/***************** S O N C h a n M a x T i m e ( ) ******************
** SONChanMaxTime:-  returns the end time contained in the last data
** block of the channel by looking in the file header
** Return 0 if channel is off or deleted
*/
SONAPI(TSTime) SONChanMaxTime(short fh, WORD chan)
{
    assert(chan < SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if ( (SONChanKind (fh, chan) == ChanOff) )
        return 0;
    else                                /* get last block */
        return(SONChanPnt(fh,chan)->maxChanTime );
}

/****************** S O N I n t l M a x T i m e ( ) *************************
** SONIntlMaxTime:-  returns the maximum time in the son file - does this
** by looking for the maximum time in any of the channels
*/
SONAPI(TSTime) SONIntlMaxTime(short fh)
{
    WORD i;
    long max = 0;
    long chanMax = 0;
    
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    for (i = 0; i < SONMAXCHANS; ++i)   /* find max time of each channel */
    {
        chanMax = SONChanMaxTime(fh, i);
        if (chanMax > max)
            max = chanMax;              /* store new max */
        else if (chanMax < 0)
            return chanMax;             /* return error if there is one */
    }

    return max;
}

/****************** S O N M a x T i m e ( ) *************************
** SONMaxTime:-  returns the maximum time in the son file - does this
** by reading from the file header
*/
SONAPI(TSTime) SONMaxTime(short fh)
{
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    return (files[fh].headP->maxFTime );
}

/****************** S O N U p d a t e M a x T i m e s ( ) *************
**
** SONUpdateMaxTimes :- Updates the stored max time information in the
** file headers to ensure it is correct.
*/
SONAPI(long) SONUpdateMaxTimes(short fh)
{
    TpFileHead  fileHeadP;                  /* pointer to file header block */
    TpChannel   channelWkP;                   /* a working ptr for channels */
    WORD        wLoop;                        /* channel counter            */

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    fileHeadP = files[fh].headP;             /* equivalenet of Pascal WITH */
    for (wLoop=0, channelWkP=files[fh].chanP;
            wLoop<SONMAXCHANS;
            wLoop++, channelWkP++)
        channelWkP->maxChanTime = SONIntlChanMaxTime(fh, wLoop);
    fileHeadP->maxFTime = SONIntlMaxTime(fh);

    return(0);     /* no errors */ 
}


/****************** S O N E x t e n d M a x T i m e ( ) *************
**
** SONExtendMaxTime :- will set the file max time to be the new value,
** provided it is greater than the current max time
*/
SONAPI(void) SONExtendMaxTime(short fh, long time)
{
    TpFileHead  fileHeadP;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    fileHeadP = files[fh].headP;
    if (time > fileHeadP->maxFTime)
    {
        fileHeadP->maxFTime = time;
        files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
    }
}


/****************** S O N G e t F i r s t D a t a ( ) ***************
**
** SONGetFirstData :- returns the offset to where the first data
** in the file is written, ie after headers.
*/
SONAPI(long) SONGetFirstData(short fh)
{
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    return files[fh].headP->firstData;
}


/*************** S O N G e t N e w F i l e N u m ( ) *********************
**
**  Search through the array of files, and return the index of the first
**  available file num (ie which doesnt correspond to an open file). Calling
**  this routine does not mark the corresponding file as allocated in anyway.
**
**  This internal routine is called by SONOpenOldFile, SONOpenNewFile to get
**  an index to an unopened file.
**
*/
SONAPI(short) SONGetNewFileNum (void)
{
    short i;

    for (i = 0; i < MAXFILES; ++i)
    {
        if (!files[i].opened)
            return i ;
    }
    return SON_OUT_OF_HANDLES;
}

/****************** S O N F i l e H a n d l e ( ) ********************
**
** SONFileHandle :- returns the file handle - to allow complete bypass
**                  of all the SON routines - if you must...
*/
#if defined(macintosh) || defined(_MAC)
SONAPI(int) SONFileHandle(short fh)
#else
    #ifdef LLIO
    SONAPI(int) SONFileHandle(short fh)
    #else
    SONAPI(FILE *) SONFileHandle(short fh)
    #endif
#endif
{
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);

#if defined(macintosh) || defined(_MAC)
    return files[fh].refNum;
#else
    #ifdef LLIO
        return files[fh].handle;
    #else
        return files[fh].handle;
    #endif
#endif
}

/****************** S O N C a n W r i t e ***************************
**
** Returns TRUE if OK to write to this file
*/
SONAPI(BOOLEAN) SONCanWrite(short fh)
{
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    return !files[fh].bReadOnly;
}

/****************** S O N I n i t F i l e s ( ) **********************
**
** SONInitFiles:-  Used to set all the files to completely undefined
** This should only be called before the files have been used, otherwise
** you will lose memory allocated to headP and to chanP.
*/
SONAPI(void) SONInitFiles(void)
{
    int i;
    for (i=0; i<MAXFILES; i++)
    {
        files[i].opened = FALSE;    /* set file not open */
        files[i].defined = FALSE;   /* nor are pointers defined */
        files[i].buffSet = FALSE;   /* no buffer yet defined */
        files[i].bReadOnly = FALSE; /* can be set TRUE on open of old file */
        files[i].lastchanRead = 0xffff; /* nothing read yet */
    }
}

/**************** S O N S e t B u f f S p a c e ( ) **********************
**
** SONSetBufSpace:- Used to set space for the biggest buffer needed for the
** file.  Search the channels for the biggest space and allocate a disk
** buffer of that size of DISKBLOCK.. whichever is the larger.
** If we fail, return SON_NO_FILE or SON_OUT_OF_MEMORY.  If the buffer is already
** allocated, we release the memory, assuming that a size change is wanted
** by the user.
*/
SONAPI(short) SONSetBuffSpace(short fh)
{
    WORD size;
    short i;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if ((fh >= MAXFILES) || (!files[fh].defined))
        return SON_NO_FILE;
    size = DISKBLOCK;                /* always space for one disk block */
    for (i=0; i<SONMAXCHANS; i++)    /* find maximum size required */
    {
        WORD phySz = SONChanPnt(fh, i)->phySz;     /* get size wanted */
        if (phySz>size)
            size = phySz;
    }
                         /* if we dont have a correct buffer already */
    if (!files[fh].buffSet || (files[fh].headP->bufferSz != size))
    {
#ifdef USEHANDLES
        if (files[fh].buffSet)
        {
            HFree(files[fh].bufferH);    /* release used memory */
            files[fh].buffSet = FALSE;
        }
                                         /* get new buffer space */
        files[fh].bufferP = (TpDataBlock ) HAlloc( size, &files[fh].bufferH);
#else
        if (files[fh].buffSet)           /* if pointer already set... */
        {
            F_free(files[fh].bufferP);   /* ...release used memory */
            files[fh].buffSet = FALSE;
    
        }                                /* and get wanted size*/
        files[fh].bufferP = (TpDataBlock ) F_malloc( size );
#endif
        if (!files[fh].bufferP)
        {
            files[fh].buffSet    = FALSE ;
            return SON_OUT_OF_MEMORY;        /* no space left, so give up */
        }
        files[fh].headP->bufferSz = size;/* remember size of buffer */
        files[fh].buffSet = TRUE;
    }
    return 0;
}

/****************** S e t F i l e S p a c e ( ) ********************
**
** Routine to set up space in memory for the file header and for the
** channel table.  This is only done if the memory has not already
** been allocated.  We do not yet set space for the file buffer as it
** is not yet known how big the buffer must be.  We return 0 if all was
** OK.  Otherwise we return SON_OUT_OF_MEMORY or SON_NO_FILE
**
** Once we have allocated space, we fill in as much as we can for the
** file header and the channels.
*/
static short SetFileSpace(short fh, WORD extra)
{
    TpFileHead fileHeadP;               /* will be pointer to the file head */
    TpChannel cP;                                    /* points at a channel */
    TSpeedPtr FAR* speedP;                /* to fill in the speed pointers */
    WORD hsize;                                /* for size of the file head */
    int i;

    assert(fh >= 0);
    assert(fh < MAXFILES);

    if (fh >= MAXFILES)     /* check a real file */
        return SON_NO_FILE;     /* error if not legal */

/* get memory space to hold the file header and the channel space */

    if (!files[fh].defined) /* do nothing if already set */
    {
        hsize = sizeof(TFileHead);  /* should be same as DISKBLOCK */
#ifdef USEHANDLES
        files[fh].headP = (TpFileHead ) HAlloc(hsize, &files[fh].headH);
        if (!files[fh].headP)
            return SON_OUT_OF_MEMORY;
            
        files[fh].chanP = (TpChannel) HAlloc(CHANSIZE, &files[fh].chanH);
        if (!files[fh].chanP)
        {
            HFree(files[fh].headH); /* release used memory */
            return SON_OUT_OF_MEMORY;   /* and say bad luck */
        }
        
        files[fh].speedP = (TpSpeedPtr)HAlloc(SONMAXCHANS * sizeof(TSpeedPtr), 
                                        &files[fh].speedH);
        if (!files[fh].speedP)
        {
            HFree(files[fh].headH); /* release used memory */
            HFree(files[fh].chanH);
            return SON_OUT_OF_MEMORY;   /* and say bad luck */
        }
#else
        files[fh].headP = (TpFileHead ) F_malloc( hsize);/* get space */
        if (!files[fh].headP)       /* returns 0 if failed */
            return SON_OUT_OF_MEMORY;
            
        files[fh].chanP = (TpChannel) F_malloc( CHANSIZE );
        if (!files[fh].chanP)       /* did we get space? */
        {
            F_free(files[fh].headP);  /* release used memory */
            return SON_OUT_OF_MEMORY;   /* and say bad luck */
        }
        
        files[fh].speedP = (TpSpeedPtr) F_malloc(SONMAXCHANS * 
                            sizeof(TSpeedPtr));
        if (!files[fh].speedP)      /* did we get space? */
        {
            F_free(files[fh].headP);  /* release used memory */
            F_free(files[fh].chanP);
            return SON_OUT_OF_MEMORY;   /* and say bad luck */
        }
#endif
        files[fh].defined = TRUE;       /* say job done */
    }

/* Now fill in as much as we can of file header */

    fileHeadP = files[fh].headP;        /* take local copy */
    fileHeadP->systemID = REVISION;     /* system version */
    fileHeadP->channels = SONMAXCHANS;  /* max number of channels in file */
    extra = (WORD)ROUND_TO_DB(extra);         /* round up to a disk block */
    fileHeadP->extraData = extra;       /* assign extra space in file */
    fileHeadP->firstData = DISKBLOCK+extra+CHANSIZE;/* points to first data */
    fileHeadP->chanSize = CHANSIZE;     /* size of the channel area */
    fileHeadP->fileState = NormalWrite; /* start with normal writing */
    F_memcpy(fileHeadP->copyright,COPYRIGHT,LENCOPYRIGHT);
    F_memcpy(fileHeadP->serialNum,"00000000",LENSERIALNUM);
    fileHeadP->osFormat = OSFORMAT;     /* say which format for byte swap */
    fileHeadP->maxFTime = 0;            /* no data in file yet */

    for (i=SON_NUMFILECOMMENTS-1; i>=0; i--)
        fileHeadP->fileComment[i].len=0;         /* set zero length strings */

                      /* now set all channels to be unused and fill in info */

    for (i=SONMAXCHANS, cP=files[fh].chanP, speedP=files[fh].speedP;
            i ; i--, cP++, speedP++)
    {
        cP->kind = ChanOff;         /* say not in use */
        cP->pad  = 0;               /* zero the rest of the word */
        cP->delSize = 0;            /* no deleted blocks yet */
        cP->nextDelBlock = CHAINEND;/* so no deleted chain */
        cP->firstBlock = CHAINEND;  /* undefined chain */
        cP->lastBlock = CHAINEND;
        cP->blocks = 0;             /* no data blocks yet */
        cP->nExtra = 0;             /* no attached data */
        cP->preTrig = 0;            /* no pretrigger */
        cP->phySz = 0;              /* no size set */
        cP->maxData = 0;            /* no number of items set */
        cP->maxChanTime = 0;        /* no data yet */
        cP->phyChan = 0;            /* no channel set */
        cP->v.adc.scale = (float)1.0; /* set default adc scale */
        cP->v.adc.offset = (float)0.0; /* seems waste of space */
        cP->idealRate = (float)100.0;
        speedP->prevBlock= CHAINEND;       /* no previous block */
        str2lstr(&cP->v.adc.units," Volt",SON_UNITSZ);
        str2lstr(&cP->comment,"No comment",SON_CHANCOMSZ); /* set comment */
        str2lstr(&cP->title,"untitled",SON_TITLESZ); /* no title */
    };
    return 0;
}


/******************** S O N G e t F r e e C h a n ( ) *********************
**
** Find an unused data channel. Either return it, or the SON_NO_CHANNEL error.
*/
SONAPI(short) SONGetFreeChan(short fh)
{
    WORD i;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    for (i=0; i<SONMAXCHANS; i++)      /* search forwards */
    if (SONChanKind(fh,i) == ChanOff)
        return (short)i;
    return SON_NO_CHANNEL;
}

/******************** S O N S e t A D C C h a n ( ) **********************
**
** fh       The handle for the file
** chan     The channel in the file to use, must not be used.
** adcChan  The physical channel from which the data came
** dvd      The divide down from the basic ADC rate in the file
**          head for this channel.
** buffSz   The physical size, in bytes, of the disk buffer to be
**          used.  This MUST be a multiple of DiskBlock.
** com      The user comment for the channel
** title    The channel title for the channel
** ideal    The ideal sampling rate for the channel (wanted by user)
** scl,offs The data *scale +offset gives real units
** unt      String holding the actual units to use.
**
** We do not change the block pointers as these have already been set to
** -1 (undefined) or the deleted chain is already set.
*/
SONAPI(short) SONSetADCChan(short fh, WORD chan, short adcChan, short dvd, short buffSz,
      TpCStr com, TpCStr title, float ideal, float scl, float offs, TpCStr unt)
{
    TpChannel cP;             /* will be pointer to the channel */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (chan>=SONMAXCHANS)
        return SON_NO_CHANNEL;
    if (files[fh].bReadOnly)
        return SON_READ_ONLY;
    cP = SONChanPnt(fh, chan);/* get a pointer to our channel */
    if (cP->kind != ChanOff)
        return SON_CHANNEL_USED;    /* must be free channel */

    cP->kind = Adc;           /* channel is now booked for use */
    cP->phySz = buffSz;       /* actual size per block on disk */
    cP->maxData = (WORD)((buffSz-SONDBHEADSZ)/sizeof(TAdc)); /* vals per block*/
    cP->phyChan = adcChan;    /* the actual input used */
    cP->idealRate = ideal;    /* Ideal sampling rate */
    cP->v.adc.scale = scl;    /* to translate to real units */
    cP->v.adc.offset = offs;
    cP->v.adc.divide = dvd;
    str2lstr(&cP->v.adc.units,unt,SON_UNITSZ);
    str2lstr(&cP->comment,com,SON_CHANCOMSZ); /* set the comment */
    str2lstr(&cP->title,title,SON_TITLESZ);   /* and the title */

    files[fh].updateHead = TRUE;    /* head needs updating now */
    files[fh].speedP[chan].prevBlock = CHAINEND;   /* no previous pointers */
    return 0;
}

/****************** S O N S e t A D C M a r k C h a n ( ) ****************
**
** fh       The handle for the file
** chan     The channel in the file to use, must not be used.
** adcChan  The physical channel from which the data came
** dvd      The divide down from the basic ADC rate in the file
**          head for this channel.
** buffSz   The physical size, in bytes, of the disk buffer to be
**          used.  This MUST be a multiple of DiskBlock.
** com      The user comment for the channel
** title    The channel title for the channel
** ideal    The ideal sampling rate for the channel (wanted by user)
** scl,offs The data *scale +offset gives real units
** unt      String holding the actual units to use.
** points   Number of data points in each item.
** preTrig  ADC points occuring before trigger point
**
** This just uses SetADCChan, then appends marker specific actions
*/
SONAPI(short) SONSetADCMarkChan(short fh, WORD chan, short adcChan, short dvd,
                 short buffSz, TpCStr com, TpCStr title, float ideal, float scl,
                            float offs, TpCStr unt, WORD points, short preTrig)
{
    short   err;
    TpChannel cP;             /* will be pointer to the channel */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    err=SONSetADCChan(fh,chan,adcChan,dvd,buffSz,com,title,ideal,scl,offs,unt);
    if (err != 0)
        return (err);
    cP = SONChanPnt(fh, chan);    /* get a pointer to our channel */
    cP->kind = AdcMark;
    cP->nExtra = (WORD)(points << 1);   /* Add AdcMark info to channel data */
    cP->preTrig = preTrig;
    cP->maxData = (WORD)((buffSz-SONDBHEADSZ) / (sizeof(TMarker)+(points<<1)));
    return 0;
}

/*************** S O N S e t R e a l M a r k C h a n ( ) **************
**
** fh       The handle for the file
** chan     The channel in the file to use, must not be used.
** phyChan  The physical channel from which the data came
** buffSz   The physical size, in bytes, of the disk buffer to be
**          used.  This MUST be a multiple of DiskBlock.
** com      The user comment for the channel
** title    The channel title for the channel
** ideal    The ideal sampling rate for the channel (wanted by user)
** min,max  The probable minimum and maximum data values
** unt      String holding the actual units to use.
** points   Number of data points in each item.
**
** This just uses SetADCChan, then appends real marker specific actions
*/
SONAPI(short) SONSetRealMarkChan(short fh, WORD chan, short phyChan,
                   short buffSz, TpCStr com, TpCStr title, float ideal,
                               float min, float max, TpCStr unt, WORD points)
{
    short   err;
    TpChannel cP;             /* will be pointer to the channel */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    err = SONSetADCChan(fh, chan, phyChan, 1, buffSz, com, title, ideal,
                                                 1.0F, 0.0F, unt);
    if (err != 0)
        return (err);
    cP = SONChanPnt(fh, chan);          /* get a pointer to our channel */
    cP->kind = RealMark;
    cP->nExtra = (WORD)(points << 2);  /* Add RealMark info to channel data */
    cP->maxData = (WORD)((buffSz-SONDBHEADSZ) / (sizeof(TMarker)+cP->nExtra));
    cP->v.real.min = min;     /* expected range of values */
    cP->v.real.max = max;

    return 0;
}

/***************** S O N S e t T e x t M a r k C h a n ( ) ****************
**
** fh       The handle for the file
** chan     The channel in the file to use, must not be used.
** phyChan  The physical channel from which the data came.
** buffSz   The physical size, in bytes, of the disk buffer to be
**          used.  This MUST be a multiple of DiskBlock.
** com      The user comment for the channel
** title    The channel title for the channel
** ideal    The ideal sampling rate for the channel (wanted by user)
** unt      String holding the actual units to use.
** points   Number of characters appended to each item (includes null).
**
** This just uses SetADCChan, then appends text marker specific actions
*/
SONAPI(short) SONSetTextMarkChan(short fh, WORD chan, short phyChan,
                                  short buffSz, TpCStr com, TpCStr title,
                                    float ideal, TpCStr unt, WORD points)
{
    short   err;
    TpChannel cP;             /* will be pointer to the channel */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    err = SONSetADCChan(fh, chan, phyChan, 1, buffSz, com, title, ideal,
                                                 (float)1.0, (float)0.0, unt);
    if (err != 0)
        return (err);
    cP = SONChanPnt(fh, chan); /* get a pointer to our channel */
    cP->kind = TextMark;
    cP->nExtra = points;       /* Add TextMark info to channel data */
    cP->maxData = (WORD)((buffSz-SONDBHEADSZ) / (sizeof(TMarker) + points));

    return 0;
}

/****************** S O N S e t I n i t L o w ( ) *****************
**
** Set the initLow datum for an eventBoth channel
*/
SONAPI(void) SONSetInitLow(short fh, WORD chan, BOOLEAN bLow)
{
    TpChannel cP;             /* will be pointer to the channel */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(SONChanKind(fh, chan) == EventBoth);

    cP = SONChanPnt(fh, chan);              /* get a pointer to our channel */
    if (cP->kind != EventBoth)                      /* Only EventBoth chans */
        return;
    if (cP->blocks != 0)  /* Cannot set initLow after data has been written */
        return;                       /* as it will not accomplish anything */
    cP->v.event.initLow = bLow;                   /* copy values, no checks */
    files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
    return;
}

/***************** S O N S e t E v e n t C h a n ( ) ****************
**
** used to associate a given channel with an event data channel. Returns
** 0 if all OK, otherwise an error code.
**
** fh           The handle for the file
** chan         The channel in the file to use, must not be inUse.
** evtChan      The physical channel from which the data came
** buffSz       The physical size, in bytes, of the disk buffer to be
**              used.  This MUST be a multiple of DiskBlock.
** com          The user comment for the channel
** title        The channel title for the channel
** ideal        The ideal sampling rate for the channel (estimated by user)
** evtKind      The actual event type used
**
** We do not change the block pointers as these have already been set to
** -1 (undefined) or the deleted chain is already set.
*/
SONAPI(short) SONSetEventChan(short fh, WORD chan, short evtChan, short buffSz,
                       TpCStr com, TpCStr title, float ideal, TDataKind evtKind)
{
    TpChannel cP;                         /* will be pointer to the channel */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (chan>=SONMAXCHANS)
        return SON_NO_CHANNEL;             /* must be legal channel number */
    if (files[fh].bReadOnly)
        return SON_READ_ONLY;
    cP = SONChanPnt(fh, chan);             /* get a pointer to our channel */
    if (cP->kind != ChanOff)
        return SON_CHANNEL_USED;           /* must be free channel */

    cP->kind = evtKind;                    /* this books the channel */
    cP->phySz = buffSz;                    /* size per block on disk */
    cP->maxData = (WORD)((buffSz-SONDBHEADSZ) / ( (evtKind==Marker) ?
              sizeof(TMarker) : sizeof(TSTime) ) );   /* per buff */
    cP->phyChan = evtChan;
    cP->idealRate = ideal;                 /* Ideal sampling rate */
    cP->v.event.initLow = FALSE;           /* Force known initial state */
    str2lstr(&cP->comment,com, SON_CHANCOMSZ); /* copy channel comment */
    str2lstr(&cP->title,title,SON_TITLESZ);    /* and the title */

    files[fh].updateHead = TRUE;           /* head needs updating now */
    files[fh].speedP[chan].prevBlock = CHAINEND; /* no previous pointers */
    return 0;
}

/**************** F i l l I n L i n k s ( ) ************************
**
** For every used channel, fill in the forward references if the file was
** written with FastWrite.
**
** If converting, should this be done???dont really know
*/
static short FillInLinks(short fh)
{
    TpChannel cP;
    int i;
    short err;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (files[fh].bReadOnly)
        return SON_READ_ONLY;
    cP = SONChanPnt(fh, 0);              /* start with first channel */
    for (i=SONMAXCHANS; i; i--)
        if (cP->kind != ChanOff)       /* only if channel is on */
        {
            long next = CHAINEND;            /* last chan link undefined */
            long here = cP->lastBlock; /* final block of channel */
            if (cP->blocks)            /* if there is some data */
            {
                do
                {
                    err = SONSetSucc(fh,here,next);/* fill in the chain */
                    if (err)
                        return err;          /* give up on error */
                    next = here;             /* prepare to move on */
                    here = workP->predBlock; /* chase backwards */
                }while (here != CHAINEND);   /* until the start */
            }
        }
    return 0 ;
}

/***************** C h e c k C o p y r i g h t ( ) ******************
**
** Given a copyright string, checks it is OK. Used by OpenOldFile.
** Returns TRUE if all OK. Replacement for strncmp, which doesn't
** exist for Windows.
*/
static BOOLEAN CheckCopyright (TpCStr copyright)
{
    int i;

    for (i=0; i<LENCOPYRIGHT; i++)
    {
        if (copyright[i] != COPYRIGHT[i])
            return (FALSE);
    }
    return (TRUE);
}

/***************** S O N O p e n O l d F i l e ( ) ***********************
**
** SONOpenOldFile:- Open up a file in the current file space.  Given a name
** we read the file header, check revision number, read header, then
** create space for channels, read channels.  Finally create space for
** the buffer.  We return 0 if OK, or an error code.
**
** The routines SONInitFiles and SONSetActiveFile must be called prior to this,
** otherwise the system will not be set up.  This routine can be called
** at any other time after this.  Any open file is shut down first.
**
** Changes that NJC would like to make - first call SONGetNewFileNum before
** calling SONSetActiveFile and SONOpenOldFile, if it fails then dont call
** SONOpenOldFile. Ie the user should explicitly close a file before a new one
** can be opened
**
** The Mac version expects a pascal string, ie a length byte followed by the
** string for the file name
**
** If we are converting a file then open it read only, as we do not
** want to allow any accidental writes to the wrong file
**
** 20.6.91 This procedure is being split up into 2 - SONOpenOldFile and
** SONReopenOldFile ReopenOldOld file expects to be called additionally with an
** open file descriptor and will not try an open the file again. This has been
** introduced for the benefit of MacApp programmers, where MacApp tries to open
** the file and  generates a file descriptor. This saves me having to override
** MacApp. For other users of the library they can still call SONOpenOldFile as
** before, and need not worry about it. But now SONOpenOldFile just does an
** open, then passes the file descriptor/ file handle to SONReopenOldFile
**
** On the Mac, if fileRefNum is <= 0, return SON_NO_FILE
** Here, if file refnum is -1 on Mac, or NULL in C or whatever, just return it
** straight back. It is up to the caller to check that the file has opened
** properly.
*/
#if defined(macintosh) || defined(_MAC)
SONAPI(short) SONReopenOldFile(short fileRefNum, BOOLEAN bReadOnly)
#endif

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
#ifdef LLIO
SONAPI(short) SONReopenOldFile(int file, BOOLEAN bReadOnly)
#else
SONAPI(short) SONReopenOldFile(FILE *file, BOOLEAN bReadOnly)
#endif
#endif

{
    short           i;          /* used as the result code */
    TpFileHead fileHeadP;       /* equivalenet of Pascal WITH */
    TpChannel cP;         /* same for the channels */
    TpChannel channelWkP;       /* a working pointer for the channels*/
    WORD            adcCount;   /* used to count ADC chans on old files */
    BOOLEAN         system2;
    BOOLEAN         system3;
    BOOLEAN         system4;
    short           fh;

#if defined(macintosh) || defined(_MAC)
    /*info we dont really need */
#if qDebug2
    HVolumeParam hRec;
    WDPBRec tempRec;
    DirInfo infRec;
    char volName[100];

    tempRec.ioCompletion = NULL;
    tempRec.ioNamePtr = volName;
    i = PBHGetVol ((WDPBPtr)&tempRec, false);
    if (i != noErr)
        return i;

    fprintf(stderr, "Default volume: %P\n", volName);
    fprintf(stderr, "vRefNum: %d\n", tempRec.ioVRefNum);
    fprintf(stderr, "vRefNum of default directory: %d\n",tempRec.ioWDVRefNum);

    hRec.ioCompletion = NULL;
    hRec.ioNamePtr = volName;
    hRec.ioVRefNum = tempRec.ioWDVRefNum;
    hRec.ioVolIndex = 0;
    i = PBHGetVInfo((HParmBlkPtr)&hRec, false);
    if (i != noErr)
        return i;
    fprintf(stderr, "Volume of default directory: %P\n", volName);

    fprintf(stderr, "dirID of default directory: %d\n", tempRec.ioWDDirID);
    infRec.ioCompletion = NULL;
    infRec.ioNamePtr = volName;
    infRec.ioVRefNum = tempRec.ioWDVRefNum;
    infRec.ioFDirIndex = -1;
    infRec.ioDrDirID = tempRec.ioWDDirID;
    i = PBGetCatInfo((CInfoPBPtr)&infRec, false);
    if (i != noErr)
        return i;
    fprintf(stderr, "Default directory: %P\n", volName);
#endif

    i = SetFPos (fileRefNum, fsFromStart, 0); /* check file really is open */
    if (i)
        return i;                  /* we probably have a bad ref num */
#endif

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
/*
** check we really have opened the file. Maybe do a seek, like on Mac to
** find out. Paul C can decide on this one
*/
#ifdef LLIO
    if (file <= 0)                 /* Can we check for silly handle */
#else
    if (file == NULL)
#endif
        return SON_NO_FILE ;           /* we have been passed a dodgy file ref */
#endif

    if ((fh=SONGetNewFileNum())<0) /* Get a file handle or fail */
        return fh ;

    i = SetFileSpace(fh,0);        /* allocate required space */
    if (i)
        return i;                  /* space not allocated */
#if defined(macintosh) || defined(_MAC)
    files[fh].refNum = fileRefNum; /* save file ref num */
#endif

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
    files[fh].handle = file;       /* save the file stream or handle */
#endif
    files[fh].bReadOnly = bReadOnly;    /* record file mode */
    fileHeadP = files[fh].headP;   /* to save repeated lookups */
    cP = files[fh].chanP;    /* ditto */
    files[fh].opened = TRUE;       /* we have a file open */
    files[fh].updateHead = FALSE;  /* no need for update, yet */

    i = SONRead(fh, fileHeadP, DISKBLOCK, 0L); /* read the file header */
    if (i != 0)                                          /* if any error... */
        goto checkClose;           /* ...report it and clean up */

    /*
    ** Now check that this is indeed likely to be one of our files
    ** otherwise we must give an error.
    */
    adcCount = 0;                  /* no adc channels yet */
    system2  = FALSE;              /* not a version 2 file */
    system3  = FALSE;              /* not a version 3 file */
    system4  = FALSE;              /* not a version 4 file */

#ifdef SONCONVERT
    i = InitConvert(fileHeadP);    /* check which way to convert */
    if (i != 0)                    /* the file, and set 'convert' variable */
        goto checkClose;
    if ((fileHeadP->osFormat != MACFORMAT) && 
        (fileHeadP->osFormat != DOSFORMAT) &&
        (fileHeadP->systemID != 3) &&
        (fileHeadP->systemID != 2))
#else
    if ((fileHeadP->osFormat != OSFORMAT) &&
        (fileHeadP->systemID != 3) && 
        (fileHeadP->systemID != 2))
#endif
    {
#if qDebug
        fprintf (stderr, "ERROR: osFormat field is wrong in hex %x\n",
            fileHeadP->osFormat);
#endif
        i = SON_WRONG_FILE;
        goto checkClose;
    }
#ifdef SONCONVERT
    gOrigVersion    = fileHeadP->systemID ; /* before we upgrade it */
#endif
    if ((fileHeadP->systemID==1) || (fileHeadP->systemID==2))
    {
        if (fileHeadP->systemID==1)
            fileHeadP->extraData = 0; /* convert sys 1 to sys 2 */
        fileHeadP->systemID = 3;   /* update header in case written back */
        system2 = TRUE;            /* will need system 2 processing to 3 */
    }

    if (fileHeadP->systemID==3)    /* now bring up to system 4 */
    {
        fileHeadP->systemID = 4;   /* update header in case written back */
        fileHeadP->osFormat = OSFORMAT; /*** convert sys 3 to sys 4 */
        system3 = TRUE;            /* will need system 2 processing to 3 */
    }

    if (fileHeadP->systemID==4)    /* now bring up to system 5 */
    {
        fileHeadP->systemID = 5;   /* update header in case written back */
        system4 = TRUE;            /* will need system 4 processing to 5 */
    }
#ifdef SONCONVERT
    sprintf (gMessageStr,"version %d, ID %d,copyright %*s", gOrigVersion,
                fileHeadP->systemID, LENCOPYRIGHT, fileHeadP->copyright);
#if qDebug
    fprintf (stderr, "%s\n",gMessageStr);
#endif
#endif

                  /* must have correct revision level and copyright message */
    if ((fileHeadP->systemID!=REVISION) || (!CheckCopyright(fileHeadP->copyright)))
    {
        i = SON_WRONG_FILE;
        goto checkClose;
    }

    /*
    ** read in the table of channels, assuming the correct size. We correct
    ** the divide down on old files. Note that old files could not generate
    ** their own adc channels, so there should not be problems with wrong
    ** counts of channels.
    */
    i = SONRead(fh, cP, CHANSIZE, DISKBLOCK);/* read channel info */
    if (i)
        goto checkClose;           /* report error and clean up */

#ifdef SONCONVERT
    if (gConvert == TONATIVE)
        for (i=SONMAXCHANS, channelWkP=cP; i ; i--, channelWkP++)
            ConvertChannelArray (channelWkP) ;
#endif

    for (i=SONMAXCHANS, channelWkP=cP; i ; i--, channelWkP++)
    {
        switch (channelWkP->kind)
        {
            case Adc:
                adcCount++;        /* increase count of channels */
            break;

            case ChanOff:
            {
                channelWkP->blocks = 0;  /* make sure no data! */
                channelWkP->firstBlock = CHAINEND;
                channelWkP->lastBlock = CHAINEND;
            }
            break;

            case EventRise:
            case EventFall:
            case EventBoth:
            case Marker:
            case AdcMark:
            case RealMark:
            case TextMark:
            break;

            default:
                i = SON_CORRUPT_FILE;
                goto checkClose;
        }/* switch */
    }
    if (system2)               /* we must change the Adc divide to system 3 */
        for (i=SONMAXCHANS, channelWkP=cP; i ; i--, channelWkP++)
            if (channelWkP->kind==Adc)
                channelWkP->v.adc.divide *= adcCount;
    if (system3)                   /* we must change the header to system 4 */
        for (i=SONMAXCHANS, channelWkP=cP; i ; i--, channelWkP++)
        {
            channelWkP->nExtra = 0; /* Generate correct data */
            channelWkP->preTrig = 0;
            channelWkP->free0 = 0;
        }
    if (system4)                   /* we must change the header to system 5 */
    {
        for (i=SONMAXCHANS, channelWkP=cP; i ; i--, channelWkP++)
        {
            if (channelWkP->comment.len > SON_CHANCOMSZ)    /* truncate comment */
               channelWkP->comment.len = SON_CHANCOMSZ;     /* if too long      */
            channelWkP->free1 = 0;                      /* init free space  */
        }

#ifdef SONCONVERT
                    /* if converting to Mac ie TONATIVE, dont update */
                    /* the headers if this file is the input */
                    /* (we defer this until we close the file */
                    /*  update the output file when we close it */
        if ((gConvert == TONATIVE) && (fh == gInputFile))
        {           /* do nothing, ie dont update the header */
        }
        else
            SONUpdateMaxTimes(fh);
#else
        SONUpdateMaxTimes(fh);
#endif

   }
    i = SONSetBuffSpace(fh);        /* allocate the buffers we need */

checkClose:                     /* here to release resources on an error */
    if (i)                      /* was there an error?*/
    {
        fileHeadP->fileState = NormalWrite; /* stop filling in of links */
        SONCloseFile(fh);       /* shut down and release allocated memory */
        return i;
    }
    else
        return fh;
}

/*
** Here is the new definition of SONOpenOldFile, which only opens the file
** then calls SONReopenOldFile
*/
#if defined(macintosh) || defined(_MAC)
/*
** Mac ROM routines expect a pascal string, as SFGetFile returns a pascal
** string if using SFGetFile, you pass it the str[63], so there should be space
** to stick a NULL on the end. It is very unlikely that a file name wil be 63
** chars in length, unless you have not been using sfgetfile, but have been
** using full pathnames
**
** Mac version: we need vRefNum, which we hope to get from stdFile SFGetFile
** also we want to make sure we have a pascal string. Also, the code needs more
** work if you want to use fsCurPerm to open the file as I don't know how to
** find out if the file was opened read only. GPS Jul/98
**
** For the pc if iOpenMode is 0, open read write, 1= read only, -1 = either
*/

SONAPI(short) SONOpenOldFile(ConstStr255Param name, short vRefNum, long dirID,
                    SignedByte perm)
{
    short i ;           /* used as the return code */
    short refNum;
    i = HOpen (vRefNum, dirID, name, perm, &refNum);
    if (i != noErr)
        return i;

    // This next call assumes that files are only opened as
    // readonly by using fsRdperm (==1). If you use fsCurPerm
    // (==0) this will not work. I don;t know of a way to get
    // the file permission back after you have opened the file.
    return SONReopenOldFile(refNum, perm == fsRdPerm);
}/* SONOpenOldFile */
#endif

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
SONAPI(short) SONOpenOldFile(TpStr name, int iOpenMode)
{
    BOOLEAN bReadOnly = iOpenMode > 0;
    #ifdef LLIO
    int      file;      /* to hold the file handle */
        #ifdef _IS_MSDOS_
    char     fname[70]; /* To get near variable holding string */
        #endif
    #else
    FILE *file;         /* file stream pointer */
    char     fname[70]; /* To get near variable holding string */
    #endif

    #ifdef _IS_WINDOWS_
    file = _lopen(name, bReadOnly ? OF_READ : OF_READWRITE);
    if (file == -1)                 /* failed to open */
    {
        if (iOpenMode < 0)          /* can we tolerate read only mode? */
        {
            file = _lopen(name, OF_READ);   /* try read only */
            bReadOnly = TRUE;
        }
        else
            return SON_READ_ONLY;
    }

    if (file == -1)
        return SON_NO_FILE;
    #else
        #ifdef LLIO
    F_strcpy(fname, name);          /* Get filename in near var */
    file = open(fname, bReadOnly ? O_BINARY|O_RDONLY : O_BINARY|O_RDWR);
    if (file == -1)                 /* failed to open */
    {
        if (iOpenMode < 0))         /* can we tolerate read only mode? */
        {
            file = open(fname, O_BINARY|O_RDONLY);   /* try read only */
            bReadOnly = TRUE;
        }
        else
            return SON_READ_ONLY;
    }

    if (file == -1)
        return SON_NO_FILE;
        #else
    F_strcpy(fname, name);          /* Get filename in near var */
    file = fopen(fname, bReadOnly ? "rb" : "r+b")
    if (!file)                      /* failed to open */
    {
        if (iOpenMode < 0))         /* can we tolerate read only mode? */
        {
            file = fopen(fname, "rb")/* try read only mode */
            bReadOnly = TRUE;
        }
        else
            return SON_READ_ONLY;
    }

    if (!file)
        return SON_NO_FILE;         /* attempt to open */
        #endif /*if LLIO else */
    #endif /* if Windows else */
    return SONReopenOldFile (file, bReadOnly);
}/* SONOpenOldFile */
#endif


/******************** S O N E m p t y F i l e ( ) ****************************
** Deletes all the data from the file (leaving channel headers the same),
** and truncates the file size, ie makes the file look like it has just been
** created. 
*/
SONAPI(short) SONEmptyFile(short fh)
{
    WORD chan;
    TpChannel cP;             /* will be pointer to the channel */
    long offset;              /* the new end of file */
    short err = 0;
    
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (files[fh].bReadOnly)
        return SON_READ_ONLY;

    for (chan = 0; chan < SONMAXCHANS; chan++)
    {
        cP = SONChanPnt(fh, chan); /* get a pointer to the channel */

        cP->delSize = 0;       /* number of blocks in deleted chain, 0=none */
        cP->nextDelBlock = CHAINEND;  /* if deleted, first block in chain */
        cP->firstBlock = CHAINEND;    /* points at first block in file */
        cP->lastBlock = CHAINEND;     /* points at last block in file */
        cP->blocks = 0;        /* number of blocks in file holding data */
        cP->maxChanTime = 0;   /* last time on this channel */

        files[fh].speedP[chan].prevBlock = CHAINEND;   /* no prev pointers */
    }
    
    files[fh].headP->maxFTime = 0;    /* max time in the data file */
    files[fh].lastchanRead = 0xffff;  /* Last chan we did a readblock from */

    offset = files[fh].headP->firstData;

#if defined(macintosh) || defined(_MAC)
    if ((err = SetEOF(files[fh].refNum, offset)) != noErr)
        return err; 
#endif

#if defined(_IS_WINDOWS_) && !defined(_MAC)
    /* Windows seek */
    if (_llseek (files[fh].handle,offset, 0) != offset)
         return SON_NO_ACCESS;
    _lwrite(files[fh].handle,NULL,0); /* truncate file to here           */
#endif  /* if Windows */

#ifdef _IS_MSDOS_
#ifdef LLIO
    {
        unsigned int wbytes;
        
        if (lseek(files[fh].handle, offset, SEEK_SET) != offset)/* LLIO read */
            return SON_NO_ACCESS;
        if (_dos_write(files[fh].handle, NULL, 0, &wbytes) != 0)
            return SON_BAD_WRITE;
    }
#else
    if (fseek (files[fh].handle, offset, SEEK_SET))     /* stdio read */
        return SON_NO_ACCESS;
    fwrite(NULL, 0, 0, files[fh].handle);
#endif  /* if LLIO else */
#endif  /* if MSDOS else */

    err = SONUpdateStart(fh);
    
    return err;
}


/***************** S O N U p d a t e S t a r t ( ) *******************
**
** This is used to write the head and channels arrays to the disk file.
** It returns 0 if Ok, else an error code.
*/
SONAPI(short) SONUpdateStart(short fh)
{
    short err;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (!files[fh].opened)
        return SON_NO_FILE;         /* must have a file */
    err=SONWrite(fh,files[fh].headP, DISKBLOCK, 0L); /* write header */
    if (err)
        return err;
    return SONWrite(fh,files[fh].chanP, CHANSIZE, DISKBLOCK);/* write chans */
}

/****************** S O N S e t F i l e C l o c k ( ) *****************
**
** Set the usPerTime and TimePerAdc values for the file
*/
SONAPI(void) SONSetFileClock(short fh, WORD usPerTime, WORD timePerADC)
{
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    files[fh].headP->usPerTime = usPerTime; /* copy values, no checks */
    files[fh].headP->timePerADC = timePerADC;
    files[fh].updateHead = TRUE;  /* Flag header to be updated on disk */
}

/***************** S O N O p e n N e w F i l e ( ) ********************
**
** SONOpenNewFile:- This is used to create a new file for writing.  To do
** this we make no assumptions about anything.  We try to create a new
** head in case it has not been set.  We also write the head to the file
** even though it will contain garbage.  We close any open handle at the
** time on this file.  fMode determines if we have fast or normal write
** to the file.
**
** CHANGES TO MAKE. This should fail if active file is already open
** It returns SON_FILE_ALREADY_OPEN, and does not close the active file
**
** On the Mac we have to explicitly create a new file before opening it
** If the file exists, we must truncate it. Then set the file type and creator
*/
#if defined(macintosh) || defined(_MAC)
SONAPI(short) SONOpenNewFile(ConstStr255Param name, short fMode, WORD extra, 
                short vRefNum, long dirID, SignedByte perm, 
                OSType creator, OSType fileType)
{
    short fileRefNum ;              /* Mac file refnum, MUST be a short */
    OSErr crErr, err ;              /* error codes from create and io ops */
    FInfo fndrInfo ;                /* for get/setting filetype and creator */
    int     i;
    TpFileHead fileHeadP;
    short   fh;

    if ((fh = SONGetNewFileNum())<0)
        return fh ;

    if ((i=SetFileSpace(fh,extra))<0)  /* try to set file header */
        return i;                   /* failed to set up header*/

/*
** Mac version is a little different, as create fails if file exists
** and open will fail if file does not exist
*/
    crErr = HCreate(vRefNum, dirID, name, creator, fileType);
    err = HOpen (vRefNum, dirID, name, perm, &fileRefNum);
    
    if (err == noErr)               /* opened the file */
    {
        if (crErr != noErr)         /* if already existed... */
        {
            if ((err = SetEOF(fileRefNum, 0)) == noErr) /* ...truncate it */
            {                                       /* get finder info */
                err = HGetFInfo(vRefNum, dirID, name, &fndrInfo); 
                if (err == noErr)
                {
                    fndrInfo.fdType = fileType;     /* Set file type */
                    fndrInfo.fdCreator  = creator ; /* and creator */
                    err = HSetFInfo(vRefNum, dirID, name, &fndrInfo); 
                }
        
            }
        }
    }
    else if (crErr != noErr)        /* check for error creating */
        err = crErr;                /* err is now our error code */
        
    if (err != noErr)
        return err ;                /* This is a Mac OSErr code*/
    files[fh].refNum = fileRefNum;  /* save the ref num if ok */
#endif

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
SONAPI(short) SONOpenNewFile(TpStr name, short fMode, WORD extra)
{
    #ifdef LLIO
    int         file;               /* local copy of file handle */
        #ifdef _IS_MSDOS_
    char        fname[70];          /* To get near variable holding string */
        #endif
    #else
    FILE        *file;              /* local copy of the file handle */
    char        fname[70];          /* To get near variable holding string */
    #endif
    short       i;
    TpFileHead  fileHeadP;
    short       fh;
    short       sErr;

    if ((fh = SONGetNewFileNum()) < 0)
        return fh ;

    if ((sErr = SetFileSpace(fh,extra)) < 0)  /* try to set file header */
        return sErr;                          /* failed to set up header*/

    #ifdef _IS_WINDOWS_
    if ( (file = _lcreat(name, 0)) < 0 )
    #else
        #ifdef LLIO
    F_strcpy(fname, name);                  /* Get filename in near var */
    if ( (file = open(fname, O_RDWR | O_CREAT | O_TRUNC | O_BINARY,
        S_IREAD | S_IWRITE)) < 0 )
        #else
    F_strcpy(fname, name);                  /* Get filename in near var */
    if (!(file = fopen(fname,"wb")))
        #endif                              /* if LLIO else */
    #endif                                  /* if Windows else */
        return SON_NO_FILE;                     /* attempt to open */
    files[fh].handle = file;                /* save the file stream */
#endif                                      /* if Msdos */

    files[fh].opened = TRUE;        /* the file is now open */
    files[fh].updateHead = TRUE;    /* the header needs writing */
    files[fh].bReadOnly = FALSE;    /* must be OK to write to it */
    fileHeadP = files[fh].headP;    /* local copy for speed */
    fileHeadP->fileState = fMode;   /* save the file write mode */

    i = SONUpdateStart(fh);           /* update the file header please */
                                    /* need to close the file if error */
    if (!i && extra)    /* Book any extra data space requested (if no err) */
        i = SONWrite(fh, workP, 1, DISKBLOCK+CHANSIZE+extra-1L);

    if (i)                          /* if errors we must close file */
    {
        fileHeadP->fileState = NormalWrite; /* stop filling in of links */
        SONCloseFile(fh);/* shut down and release any allocated memory */
        return i;                   /* return any error */
    }
    else
        return fh;                  /* return file handle */
}

/**************** S O N C l o s e F i l e ( ) ************************
**
** Close down the currently active file and fill in the backward file
** links if FastWrite was used. return 0 if all OK or an error code.
** On the Mac it is necessary to call FlushVol to write out the VolCtlBlk
** after closing. This should always (?) be done after Close on Mac but we
** need the volume ref num
*/
SONAPI(short) SONCloseFile(short fh)
{
    short       error = 0;
    int         clError = 0;                          /* mac FSclose result */
    short       vRefNum = 0;    /* only used by Mac routines, MUST BE SHORT */
    TpFileHead  fileHeadP;
    TpChannel   channelP;
    short       lowest;                     /* lowest version no we try for */
    short       i;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    fileHeadP = files[fh].headP;
    channelP  = files[fh].chanP;
    if (!files[fh].opened)
        return SON_NO_FILE; /* check open */

/*
** now see if we can downgrade the version number, we only try down to level
** 3 as its too different before that. Versions 4 and 5 are just extensions
** of version 3 with extra t ype.
*/
    lowest = 3;                             /* lowest version no we try for */
    for (i=0; i<SONMAXCHANS; i++)           /* Search through the channels  */
    {
        switch (channelP[i].kind)
        {
            case AdcMark :                    /* AdcMark requires version 4 */
            {
                if (lowest < 4)
                    lowest = 4;
                break;
            }
            case TextMark :                  /* TextMark requires version 5 */
            case RealMark :                             /* as does RealMark */
            {
                if (lowest < 5)
                    lowest = 5;
                break;
            }
            default : ;                     /* No action by default */
        }
    }

#ifdef SONCONVERT
                    /* do the updating of max time, but only */
                    /* if we didnt do it before and we had to upgrade */
                    /* from 4 to 5 */
    if ((fh == gOutputFile) && (gConvert == TONATIVE) &&
        (gOrigVersion <= 4) && (lowest >=5))
    {
        SONUpdateMaxTimes (fh) ;
        files[fh].updateHead=TRUE;          /* as head changed */
    }
#endif

    if (fileHeadP->systemID != lowest)      /* Update if we are making a change */
    {

        fileHeadP->systemID = lowest;
#ifdef SONCONVERT
                                            /* if we are converting to dos, */
                                            /* byte swap version */
    if (gConvert != TONATIVE)
        fileHeadP->systemID = lowest * 256 ;
#endif
    }

    if ((fileHeadP->fileState==FastWrite) && /* must fill in links? */
        !files[fh].bReadOnly)               /* it should never be readonly */
    {
        error = FillInLinks(fh);            /* fill in file links */
        if (error != 0)
            return error;                   /* error if failed */
        fileHeadP->fileState=NormalWrite;   /* as links filled in */

        SONUpdateMaxTimes(fh);              /* not filled in as we went */
        files[fh].updateHead=TRUE;          /* as head changed */
    }

    if (files[fh].buffSet)
    {
        files[fh].buffSet = FALSE;          /* say buffer released... */
#ifdef USEHANDLES
        HFree(files[fh].bufferH);           /* ...and release it */
#else
        F_free(files[fh].bufferP);          /* release buffer */
#endif
    }

    if (files[fh].updateHead &&             /* head needs update, and... */
        !files[fh].bReadOnly)               /* ...not read only */
        error = SONUpdateStart(fh);         /* save head on disk */

#if defined(macintosh) || defined(_MAC)
            /* if we get an error, report the first one, ie SONUpdateStart */
            /* error. But we should still try and close */
    clError = GetVRefNum (files[fh].refNum, &vRefNum);/* get volume ref Num */
    clError = FSClose (files[fh].refNum);

#if qDebug
    if (clError)
        fprintf(stderr, "ERROR: %d SONCloseFile\n", clError) ;
#endif
    FlushVol (NULL, vRefNum);               /* flush the volume anyway */
    if (!error)                             /* only return clError if*/
        error = clError ;                   /* no error earlier */

#endif

#if defined(_IS_WINDOWS_) && !defined(_MAC)
    _lclose( files[fh].handle );            /* shut the file */
#endif

#ifdef _IS_MSDOS_
#ifdef LLIO
    close( files[fh].handle );              /* shut the file */
#else
    fclose(files[fh].handle);               /* shut the file */
#endif
#endif

#ifdef USEHANDLES
    HFree(files[fh].chanH);                  /* free channel area */
    HFree(files[fh].headH);                  /* and header area */
    HFree(files[fh].speedH);                 /* and speedPtrs area */
#else
    F_free(files[fh].chanP);                   /* free the channel area */
    F_free(files[fh].headP);                   /* and the header area */
    F_free(files[fh].speedP);                  /* and speedPtrs area */
#endif
    files[fh].defined = FALSE;               /* no head space left */
    files[fh].opened = FALSE;                /* file not open */
    return error;
}

/******* S O N G e t S u c c ( )   &   S O N S e t S u c c ( ) **********
**
** SONGetSucc & SONSetSucc:- Used to get, and set the succBlock field of data
** blocks in the file.  This allows us to patch the file to string blocks
** together.  Note we use the work pointer to a single disk block.  Also
** SONGetPred, used for searches backwards through file. All of these routines
** use SONGetBlock to read a single block of data.
**
** While SON_NO_FILE == -1, SONGetSucc could return -1 for end of chain, or
** because of SON_NO_FILE error condition. 4.7.91, as yet unresolved
**
** 4.7.91 SONGetBlock no longer static, as it is used to get maxtime
*/
SONAPI(short) SONGetBlock(short fh, long offset)
{
    short err;
    
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (!files[fh].opened) return SON_NO_FILE;
    err = SONRead(fh, workP, DISKBLOCK, offset);
#ifdef SONCONVERT
        /* convert to native format if it is foreign */
    if ((gConvert == TONATIVE && fh == gInputFile) ||
        (gConvert == TOFOREIGN && fh == gOutputFile))        
        ConvertBlockHead(workP);
#endif
    return err;
}

/**************** S O N G e t S u c c ( ) ***********************
** Get the file position of the next block in the linked list
**
** fh       The handle for the file
** offset   File offset to the current block
** BUG HERE see SONGetBlock, which can return -1, so can workP->succBlock
*/
SONAPI(long) SONGetSucc(short fh, long offset)
{
    int flag;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    flag = SONGetBlock(fh,offset);  /* read block and check error */
    return (flag != 0) ? (long)flag : workP->succBlock;
}

/*************** S O N S e t S u c c ( ) **************************
** Attempt to fill in the pointer to the next block
**
** fh       The handle for the file
** offset   Points to the header of the current block
** succOffs Points to the following block
*/
SONAPI(short) SONSetSucc(short fh, long offset, long succOffs)
{
    short err;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (offset == CHAINEND)
        return 0;
    err = SONGetBlock(fh,offset);   /* get current block */
    if (err != 0)                   /* get current block */
        return err;                 /* return it if an error */
    workP->succBlock = succOffs;    /* fill in successor */
#ifdef SONCONVERT
        /* convert to foreign format if we should */
    if (gConvert == TOFOREIGN && fh == gOutputFile)      
        ConvertBlockHead(workP);
#endif
    err = SONWrite(fh,workP, DISKBLOCK, offset); /* and write it back */
#ifdef SONCONVERT
        /* convert it back again just in case it is used ... */
        /* (as it is in FillInLinks) */
    if (gConvert == TOFOREIGN && fh == gOutputFile)      
        ConvertBlockHead(workP);
#endif
    if ( err != 0)                  /* return any error */
        return err;
    files[fh].updateHead = TRUE;
    return 0;
}

/**************** S O N G e t P r e d ( ) **************************
** Get the offset to the previous block in the chain, -1 means there
** are no previous blocks.
**
** Normally only used internally.
**
** fh       The handle for the file
** offset   Offset to the current block
*/
SONAPI(long) SONGetPred(short fh, long offset)
{
    int flag;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    flag = SONGetBlock(fh,offset);
    if (flag != 0)
        return flag;               /* get requested block */
/*
** here we could do a check on the integrity of the file - the offset
** should be within the range 0 to EOF
*/
    return workP->predBlock;       /* return predecessor */
}


/**************** S O N R e a d B l o c k ( ) **************************
**
** Used to fill *workP array with data from a given position in the file
**
** fh       The handle for the file
** chan     The data channel to read from
** position The file offset to read from
*/
SONAPI(short) SONReadBlock(short fh, WORD chan, long position)
{
    short       err;
    TpSpeedPtr  pSpeed;
    TpDataBlock bufP;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (!files[fh].opened)
        return SON_NO_FILE;                                 /* must have a file */

    bufP = files[fh].bufferP;                        /* xfer buffer address */
    pSpeed = &files[fh].speedP[chan];        /* Pointer into speedup tables */
    err = 0 ;                   /* No error for now, read only if necessary */
    if ((position != pSpeed->prevBlock) ||
        (files[fh].lastchanRead != chan))
        err = SONRead(fh,files[fh].bufferP, SONChanPnt(fh, chan)->phySz, position) ;

    if (err != 0)                 /* this will preserve the return code, if */
        return err;                                  /*  it was a maccy one */

#ifdef SONCONVERT
/*
** for TONATIVE or TOFOREIGN convert the data in the block
** It only converts the block header if TONATIVE, ie if its foreign format
**  on disk
*/
    if (gConvert)
        err = ConvertDataBlock (bufP, SONChanPnt(fh, chan));
    if (err)
        return err ;
#endif

/*
** Now copy information to the speed pointers to save time in searches
*/
    pSpeed->prevBlock = position;
    pSpeed->prevStart = bufP->startTime;
    pSpeed->prevEnd   = bufP->endTime;
    pSpeed->prevPred  = bufP->predBlock;
    pSpeed->prevSucc  = bufP->succBlock;
    files[fh].lastchanRead = chan;

    return 0;
}


/****************** S O N W r i t e B l o c k ( ) ********************
**
** SONWriteBlock:- Used to write the current block defined in bufferP^ to the
** file.  We assume (so FAR) that the data is to be appended to the end of
** the file.  It would be possible to add such items in the middle, but
** this becomes very difficult indeed.
**
** We must sort out all the disk pointers, and allow for re-using deleted
** file chains.
**
** In Mac SonConvert there are 2 files open, input and output, sharing the
** same bufferP. outChannelP (only used in Mac SonConvert) is set up by a
** calling routine and is used to access the input file channel pointer
*/
SONAPI(short) SONWriteBlock(short fh)
{
    short       err ;                       /* just for Mac , unused in DOS */
    long        next;
    TpFileHead  fileHeadP;                  /* Pointer to file header */
    TpDataBlock bufP;                           /* to data buffer */
    WORD        chan;
    TpChannel   cP;
    TSTime      blkEndTime;

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (files[fh].bReadOnly)
        return SON_READ_ONLY;

    fileHeadP = files[fh].headP;
    bufP = files[fh].bufferP;
    chan = (WORD)((bufP->chanNumber & 255) - 1); /* in case event level included */
    cP = SONChanPnt(fh, chan);
    blkEndTime = bufP->endTime;
    if (cP->delSize != 0)                  /* using deleted blocks? */
    {
        next = cP->nextDelBlock;      /* next deleted block */
        cP->nextDelBlock = SONGetSucc(fh,next);     /* for next time */
        cP->delSize--;                 /* reduce blocks in the chain */
    }
    else
    {

#if defined(macintosh) || defined(_MAC)
        err = SetFPos (files[fh].refNum, fsFromLEOF, 0);     /* seek to eof */
        if (err == noErr)
            err = GetFPos (files[fh].refNum, &next);         /* get the eof */
        if (err != noErr)
            return err ;                          /* this is mac oserr code */
#endif

#if defined(_IS_WINDOWS_) && !defined(_MAC)
        next = _llseek(files[fh].handle,0L,SEEK_END);/* get end of file pos */
#endif

#ifdef _IS_MSDOS_
#ifdef LLIO
        next = lseek( files[fh].handle, 0L, SEEK_END);/* get end of file pos */
#else
        fseek(files[fh].handle,0L,SEEK_END);/* move to end of file */
        next = ftell(files[fh].handle);     /* get position */
#endif
#endif
    }

    if ((fileHeadP->fileState != FastWrite) /* if not fast write mode...*/
         && (cP->blocks != 0))              /* ...and got blocks */       
    {
        err = SONSetSucc(fh, cP->lastBlock, next);  /* fill in the link */
        if (err)
            return err;
    }

    bufP->succBlock = CHAINEND;         /* set forward and backward links */
    bufP->predBlock = (cP->blocks != 0) ? cP->lastBlock : CHAINEND;

#ifdef SONCONVERT
    if (gConvert == TOFOREIGN)
        ConvertBlockHead(bufP);
#endif
    files[fh].lastchanRead = 0xffff;    /* Stored data not valid .. */
    err=SONWrite(fh, bufP, cP->phySz, next);    /* write the data to file */

#ifdef SONCONVERT
    if (gConvert == TOFOREIGN)           /* make sure head is correct */
        ConvertBlockHead(bufP);
#endif
    if (err)
        return err;

/*
** now adjust pointers in channel area
*/
    if (cP->blocks == 0)                 /* first block? */
        cP->firstBlock = next;           /* fill in first link */
    cP->lastBlock = next;                /* fill in last block */
    cP->blocks++;                        /* increase block count */

    if (blkEndTime > cP->maxChanTime)
      cP->maxChanTime = blkEndTime;         /* update channel max time */
    if (blkEndTime > fileHeadP->maxFTime)
      fileHeadP->maxFTime =  blkEndTime;             /* and file max too */

    files[fh].updateHead = TRUE;    /* header has changed */
    return 0;                       /* done without error */
}

/************ S O N W r i t e E v e n t B l o c k ( ) *******************
**
** SONWriteEventBlock:- Output a block (or blocks) of data. If data output
** would span more than one data block it is split up into several blocks.
** This can lead to blocks of different lengths being stored on disk.
**
** fh           The handle for the file
** chan         The channel number to write this block to
** buffer       The buffer holding the data
** count        The number of data items to write
*/
SONAPI(short) SONWriteEventBlock(short fh,WORD chan,TpSTime buffer,long count)
{
    long            dataPoint = 0;     /* index of item to transfer next */
    TpChannel       cP;                  /* local channel pointer */
    TpDataBlock     bufP;                         /* local data pointer */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    cP = SONChanPnt(fh, chan);                  /* get channel pointer */
    bufP = files[fh].bufferP;                /* and data pointer */
    if (cP->kind==ChanOff)
        return SON_CHANNEL_UNUSED;

    while (count)   /* keep writing buffers until all data written */
    {
        short err;
        WORD transfer = (WORD)((count > (long)cP->maxData) ? cP->maxData : count);

        bufP->items = transfer;     /* set items in the buffer itself */
/*
** sort out which way up we are for EventBoth - I checked this, is OK. TDB
*/
        if (cP->blocks == 0)        /* if first block copy initLow */
            cP->v.event.nextLow = cP->v.event.initLow;
        bufP->chanNumber = (WORD)((cP->v.event.nextLow) ? chan+257 : chan+1);

        bufP->startTime = buffer[dataPoint];    /* first saved time */
        F_memcpy(bufP->data.int4Data, buffer+dataPoint, 
                    transfer*sizeof(TSTime));
        dataPoint += transfer;       /* move pointer on in the buffer */
        bufP->endTime = buffer[dataPoint-1];    /* last time in the buffer */
        if ((err = SONWriteBlock(fh)) != 0)     /* write block and... */
            return err;                         /* ...return any error */
        if (transfer & 1)       /* if an odd number of points invert level */
            cP->v.event.nextLow = (BOOLEAN)!cP->v.event.nextLow; 
        count -= transfer;      /* reduce the count left */
    }

    return 0;
}

/*********** S O N W r i t e M a r k B l o c k ( ) *****************
**
** SONWriteMarkBlock:- Output a block (or blocks) of data.  If data output
** would span more than one data block it is split up into several blocks.
** This can lead to blocks of different lengths being stored on disk.
**
** fh           The handle for the file
** chan         The channel number to write this block to
** buffer       The buffer holding the data
** count        The number of data items to write
*/
SONAPI(short) SONWriteMarkBlock(short fh,WORD chan,TpMarker buffer,long count)
{
    TpDataBlock bufP;   /* local data pointer */
    TpChannel cP;    /* local channel pointer */
    long dataPoint = 0;                     /* next transfer index */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    bufP = files[fh].bufferP;
    cP = SONChanPnt(fh, chan);
    if (cP->kind != Marker)
        return SON_NO_CHANNEL;

    while (count)   /* keep writing buffers until all data written */
    {
        short err;
        WORD transfer  = (WORD)((count > (long)cP->maxData) ? cP->maxData : count);

        bufP->items = transfer;     /* set items in this buffer */
        bufP->chanNumber = (WORD)(chan + 1);
        bufP->startTime = buffer[dataPoint].mark;   /* first saved time */
        F_memcpy(bufP->data.markData,buffer+dataPoint,transfer*sizeof(TMarker));/* fast copy */
        dataPoint += transfer;      /* move pointer on in buffer */
        bufP->endTime = buffer[dataPoint-1].mark;   /* last time in buffer */
        if (( err = SONWriteBlock(fh)) != 0)        /* write block and... */
            return err;                             /* ...return any error */
        count -= transfer;          /* reduce the count left */
    }

    return 0;
}

/************** S O N W r i t e E x t M a r k B l o c k ( ) **************
**
** SONWriteExtMarkBlock:- Output a block (or blocks) of data.  If the data
** would span more than one data block it is split up into several blocks.
** This can lead to blocks of different lengths being stored on disk.
** The data can be any type of extended marker - AdcMark, RealMark, TextMark.
**
** fh     The handle for the file
** chan   The channel number to write this block to
** buffer The buffer holding the data
** count  The number of data items to write
*/
SONAPI(short) SONWriteExtMarkBlock(short fh, WORD chan, TpMarker buffer, long count)
{
    TpDataBlock bufP;                                 /* local data pointer */
    TpChannel cP;                              /* local channel pointer */
    WORD nBytes;                          /* size of the data item to write */
    long dataPoint = 0;                       /* index into the data buffer */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    bufP = files[fh].bufferP;           /* Get local pointers */
    cP = SONChanPnt(fh, chan);
    if ((cP->kind != AdcMark) &&
       (cP->kind != RealMark) &&
       (cP->kind != TextMark))
        return SON_NO_CHANNEL;

    nBytes = SONItemSize(fh,chan);

#ifdef APPENDWRITE                         /* Tim's code for appending data */
    if ((cP->blocks > 0) && (ptsInLast < cP->maxItems))
    {
        short err;         /* If we do this, move these to function level ? */
        TpMarker markP;
        WORD space = (cP->maxItems - ptsInLast);
        WORD transfer = (WORD)((count > space) ? space : count);

        err = SONReadBlock(fh, chan, cP->lastBlock); /* Get back last block */
        if (err < 0)
            return err;
        if (bufP->items != ptsInLast)       /* This should never occur ... */
        {                       /* but I put it in to cause discussion */
            return SON_BAD_READ;        /* special error code ? */
        }
        markP = movePtr(bufP->data, ptsInLast*nBytes);
        bufP->items += transfer;
        F_memcpy(markP, buffer, transfer*nBytes);              /* fast copy */

        markP = movePtr(buffer, (transfer-1)*nBytes);  /* Last time pointer */
        bufP->endTime = markP->mark;                 /* last time in buffer */
        if ( (err = SONWriteBlock(fh)) != 0)          /* write block and... */
            return err;                              /* ...return any error */
        count -= transfer;                         /* reduce the count left */
    }
#endif

    while (count)            /* keep writing buffers until all data written */
    {
        short err;
        TpMarker markP = (TpMarker)movePtr(buffer, dataPoint*nBytes);   /* start addr */
        WORD transfer = (WORD)((count > (long)cP->maxData) ? cP->maxData : count);

        bufP->items = transfer;                  /* set items in the buffer */
        bufP->chanNumber = (WORD)(chan+1);   /* NB chans on disk start at 1 */
        bufP->startTime = markP->mark;                  /* first saved time */
        F_memcpy(bufP->data.markData, markP, transfer*nBytes); /* fast copy */

        dataPoint += transfer;                 /* move pointer on in buffer */
        markP = (TpMarker)movePtr(buffer, (dataPoint-1)*nBytes);       /* Last marker */
        bufP->endTime = markP->mark;                 /* last time in buffer */
        if ( (err = SONWriteBlock(fh)) != 0)          /* write block and... */
            return err;                              /* ...return any error */
        count -= transfer;                         /* reduce the count left */
    }

    return 0;
}

/**************** S O N W r i t e A d c B l o c k ( ) ******************
**
** SONWriteADCBlock:- Output a block (or blocks) of data. If the data output
** would span more than one data block it is split up into several blocks.
** This can lead to blocks of different lengths being stored on disk.
**
** fh     The handle for the file
** chan   The channel number to write this block to
** buffer The buffer holding the data
** count  The number of data items to write
** sTime  The time of the first item in the block
**
** returns the time of the next sample or -ve error code
*/
SONAPI(TSTime) SONWriteADCBlock(short fh, WORD chan, TpAdc buffer,long count,
                                                                TSTime sTime)
{
    long dataPoint = 0;     /* points at item in the buffer */
    TSTime adcTime;
    TpChannel cP;                            /* local channel pointer */
    TpDataBlock bufP;                            /* local data pointer */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    cP = SONChanPnt(fh, chan);                   /* get local pointers */
    bufP = files[fh].bufferP;
    if (cP->kind != Adc)
        return SON_NO_CHANNEL;
    adcTime = SONChanDivide(fh,chan); /* get time per adc point */

    while (count)   /* write buffers until all data written */
    {
        int err;
        WORD transfer = (WORD)((count > (long)cP->maxData) ? cP->maxData : count);
        bufP->items = transfer;         /* set items in the buffer */
        bufP->chanNumber = (WORD)(chan+1);
        bufP->startTime = sTime;        /* first saved time */
        sTime += transfer*adcTime;      /* time of next block */
        bufP->endTime = sTime-adcTime;  /* last time in block */
                         /* copy data (quickly) into the buffer for writing */
        F_memcpy(bufP->data.int2Data, buffer+dataPoint, transfer*sizeof(short));
        dataPoint += transfer;          /* move pointer on in buffer */
        if ((err=SONWriteBlock(fh))!=0) /* write block and... */
            return err;                 /* ...return any error */
        count -= transfer;              /* reduce the count left */
    }

    return sTime;
}


/**************** S O N C o m m i t F i l e ( ) ***********************
**
** SONCommitFile:- Commits a file to disk, ie ensures that file secure
** on disk. The precise actions will vary with the OS, needless to say,
** but the DOS mechanism is to update the file header on disk, followed
** by duplicating the file handle and closing it, which updates the
** directory information.
**
** fh     The handle for the file
**
** returns zero or -ve error code
*/
SONAPI(short) SONCommitFile(short fh)
{
    short err;
#if defined(macintosh) || defined(_MAC)
    short vRefNum;                      /* volume ref num */
#endif
#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
    int hand;
#endif

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    err = SONUpdateStart(fh);           /* First update the header */
    if (err < 0)
        return err;

#if defined(macintosh) || defined(_MAC)
    err = GetVRefNum (files[fh].refNum, &vRefNum);
    if (!err)
        err = FlushVol (NULL, vRefNum); /* flush data and directory info */
                                        /* to disk */
#endif

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
#ifdef LLIO                            /* DOS/Windows LLIO version          */
    hand = _dup(files[fh].handle);     /* Duplicate the handle              */
    if (hand < 0)
        return SON_NO_HANDLES;
  #ifdef _IS_WINDOWS_
    return (short)_lclose( hand );     /* shut the file for Windows         */
  #endif

  #ifdef _IS_MSDOS_
    return (short)close( hand );       /* shut the file                     */
  #endif

#else                                  /* DOS STDIO version                 */
    hand = _dup(_fileno(files[fh].handle));         /* Duplicate the handle */
    if (hand < 0)
        return SON_NO_HANDLES;
    return (short)fclose(hand);                            /* shut the file */
#endif
#endif

}


/*************** S O N F i n d B l o c k ( ) *****************************
**
** SONFindBlock:- This is used to find the first block in a file which has
** data which fits into a specified time range.  If the data overlaps
** the time range at all, it counts. -ve return is error, 0 indicates no
** data in specified range and +ve result is block offset in file.
**
** fh       The handle for the file
** chan     The data channel to search on disk
** sTime    The first time in the range we seek
** eTime    The last time in the range
**
** We return either the block offset in which the time was found, 0 if
** no time was found, or -ve for an error.  If the channel is not used
** we return SON_CHANNEL_UNUSED and if it's closed we return SON_NO_FILE.
*/
SONAPI(long) SONFindBlock(short fh, WORD chan, TSTime sTime, TSTime eTime)
{
    long block,next;            /* used as file offsets */
    BOOLEAN goForwards;         /* the search direction */
    TpSpeedPtr speedP;                               /* local pointer */
    TpChannel cP;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (!files[fh].opened)      /* is there a file? */
        return SON_NO_FILE;
    speedP = &files[fh].speedP[chan];             /* get speed pointer */
    cP = SONChanPnt(fh, chan);
    if (files[fh].headP->fileState == FastWrite)    /* useable? */
        return SON_NO_ACCESS;
    if (cP->kind==ChanOff)          /* is this channel used? */
        return SON_CHANNEL_UNUSED;
    if ((cP->blocks == 0) ||        /* any data blocks? */
                (eTime < sTime))    /* check not silly time range */
        return 0;
    /*
    ** Now we try to optimise the search using the data in speedPtrs. We
    ** begin by setting default search parameters, then try to use the
    ** info to speed up the search.
    */
    goForwards = TRUE;          /* start with forward search */
    block = cP->firstBlock;     /* set the start point */

    if (speedP->prevBlock >= 0) /* any previous reads to work on? */
    {
        if (sTime>speedP->prevEnd)  /* is wanted time after previous read? */
        {
            block = speedP->prevSucc;       /* try next block */
            if ((block<0) ||                /* end of chain, or... */
                (block > cP->lastBlock))    /* not saved yet */
            {
                block = cP->lastBlock;      /*...search from the end...*/
                goForwards = FALSE;         /*...backwards */
            }
        }
        else    /* choose from start, or last read */
        {
            if (sTime>=speedP->prevStart) /* ...in last block */
                return speedP->prevBlock; /* then we have found it */
            else
                if ((speedP->prevBlock==cP->firstBlock) &&
                        (eTime<speedP->prevStart)) /* before first */
                    return 0;

            if (sTime>speedP->prevStart / 2) /* based on time */
            {
                if (speedP->prevPred!=CHAINEND) /* somewhere to go? */
                {
                    block = speedP->prevPred;   /* search from here */
                    goForwards = FALSE;         /* backwards */
                }
                else
                    return speedP->prevBlock;   /* must exist if > sTime */
            }
        }
    }/* if previous reads to work on */

    if (goForwards) /* code to search forwards */
    {
        while (block >= 0)    /* test here as prev may adjoin fail end */
        {
            next = SONGetSucc(fh,block);    /* get 1 block into workP */
            if (workP->endTime >= sTime)    /* data in the block? */
            {
                if (workP->startTime <= eTime)
                    return block;           /* we have found the block */
                else
                    return 0;               /* block does not exist */
            }
            if (next > cP->lastBlock)       /* Past end of written data? */
                return 0;
            block = next;                   /* otherwise, on to the next */
        }
    }
    else            /* code to search backwards */
    {
        while (block >= 0)
        {                                   /* start block must be written */
            next = SONGetPred(fh,block);    /* get 1 block back into workP */

            if (workP->startTime <= sTime)  /* early enough? */
            {
                if (workP->endTime >= sTime)/* have we gone past? */
                    return block;
                else
                {                           /* not gone past */
                    block = workP->succBlock;   /* find next */
                    if ((block<0) ||            /* end of chain, or ...*/
                       (block > cP->lastBlock)) /*... not yet written */
                        return 0;               /* give up */
                    SONGetSucc(fh,block);   /* move along to next block */
                    if (workP->startTime <= eTime)
                        return block;
                    else
                        return 0;
                }
            }
            block = next;                   /* chase backwards */
        }

        block   = cP->firstBlock ;
        next    = SONGetSucc (fh,block);
        if (workP->startTime <= eTime)
            return block ;
        else
            return 0;

    }/* if goForwards */

    return (block!=CHAINEND) ? block : 0;   /* block #, error or not found */
}


/******************* S O N G e t E v e n t D a t a ( ) ******************
**
** SONGetEventData:- Used to collect event data into a buffer between two
** specific times.  The buffer is of type ADSEvents, pointing at an array
** of events.  The caller requests a number of events, the routine gives
** back the number it actually transferred, or -ve error code.
**
** fh           The handle for the file
** chan         The data channel to get times from
** evtDataP     Where to store the times we find
** max          The maximum number of events to collect
** sTime,eTime  The start and end of the period we want data in
** levLow       returned TRUE if first event in block is high going, ie level
**                  is low before first event
**
*/
SONAPI(long) SONGetEventData(short fh, WORD chan, TpSTime evtDataP, long max,
        TSTime sTime, TSTime eTime, TpBOOL levLowP, TpFilterMask pFiltMask)
{
    long block;
    long nBytes;                        /* Number of bytes in an item */
    long count = 0;                     /* number of items returned */
    BOOLEAN markerKind;
    TpSTime evtWrP = evtDataP;          /* writing pointer. */
    TpDataBlock bufP;                   /* local copy of buffer pointer */
    TpChannel cP;                       /* pointer to channel data */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(max > 0);
    assert(levLowP != NULL);
    
    bufP = files[fh].bufferP;           /* get copy of buffer pointer */
    cP = SONChanPnt(fh, chan);
    switch (SONChanKind(fh,chan))  /* check the type of the channel*/
    {
        case Marker:
        case AdcMark:
        case RealMark:
        case TextMark:
            markerKind = TRUE;
            break;

        case EventFall:
        case EventRise:
        case EventBoth:
            markerKind = FALSE;
            break;

        default:
            return SON_NO_CHANNEL;
    }

    nBytes = SONItemSize(fh,chan);           /* Size of data item, in bytes */
    if ((block = SONFindBlock(fh, chan, sTime, eTime))<=0)  /* find a block */
    {
        long err = 0;
/*
** No data in time range, but if event data have to set levLow to show state
** correctly , so look for first block after the time range requested.
** The value returned is for the first transition AFTER time range used.
*/
        if (!markerKind && (block == 0))
        {
            block = SONFindBlock(fh, chan, sTime, LONG_MAX);
            if (block > 0)                           /* If we got something */
            {
                err = SONReadBlock(fh,chan,block);
                if (err == 0)      /* level is same as first event in block */
                    *levLowP = (BOOLEAN)((bufP->chanNumber & 256) != 0);
            }
            else /* No data found, next transition WOULD be inverse of last */
            {
                block = cP->lastBlock;            /* get the last block */
                if (block > 0)
                {
                    err = SONReadBlock(fh,chan,block);
                    if (err == 0)  /* ext trans is same as after last event */
                        *levLowP = (BOOLEAN)(((bufP->chanNumber & 256) != 0) ^
                                     (bufP->items & 1));
                }
                else                                    /* no blocks at all */
                    *levLowP = cP->v.event.initLow;
            }
        }
        else
            err = block;

        return err;             /* return with error */
    }

/*
** Search through for first useful time (there must be one). Use
** different loop for markers. Work out the level for eventBoth
** data (ignored by other data types.
*/
    do
    {
        short err = SONReadBlock(fh,chan,block);/* get next data block */
        WORD dataLeft = bufP->items;            /* number of items in buff */
        if (err<0)                              /* check for read error */
            return err;

        if (markerKind)
        {
            TpMarker markRdP = bufP->data.markData; /* marker address */
            if ((count == 0) && dataLeft)     /* is this the first buffer? */
                while (markRdP->mark < sTime) /* skip unwanted data */
                {
                    markRdP = (TpMarker)movePtr(markRdP, nBytes);
                    dataLeft--;
                }

            while (dataLeft && (count<max))   /* loop to collect data */
            {
                TSTime time4 = markRdP->mark; /* get the marker time */
                if (time4 > eTime)
                    return count;              /* count of events */
                if ( !pFiltMask ||             /* is filter defined? */
                     SONFilter(markRdP, pFiltMask)) /* do we want it? */
                {
                    *evtWrP++ = time4;        /* add to output */
                    count++;
                }
                markRdP=(TpMarker)movePtr(markRdP,nBytes); /* move on pointer and */
                dataLeft--;                   /* count of data left */
            }
        }
        else
        {
            TpSTime evtRdP = bufP->data.int4Data; /* data pointer for event */
            if ((count == 0) && dataLeft)      /* is this the first buffer? */
            {
                while (*evtRdP < sTime)               /* skip unwanted data */
                {
                    evtRdP++;
                    dataLeft--;
                }
                                              /* first level */
                *levLowP = (BOOLEAN)(((bufP->chanNumber & 256) != 0) ^
                                     ((bufP->items - dataLeft) & 1));
            }

            while (dataLeft && (count < max))
            {
                TSTime time4 = *evtRdP++;
                if (time4 > eTime)
                    return count;             /* count of events */
                *evtWrP++ = time4;            /* add to output */
                count++;                      /* increase count */
                dataLeft--;                   /* decrease left in buf */
            }
        }
    /*
    ** To get here either we have filled the output (count=max) or
    ** we have emptied a buffer (dataLeft==0).
    */
        block = bufP->succBlock;      /* next block to read */
        if (block > cP->lastBlock)    /* for GRB family */
            block = CHAINEND;

    }while ( (count != max) &&        /* buffer not full yet */
             (block > 0) &&           /* and not reached eof */
             (bufP->endTime <= eTime)); /* and not off the end */

    return count;                     /* stop now! */
}

/***************** S O N G e t M a r k D a t a ( ) ***********************
** SONGetMarkData:- Used to collect marker data into a buffer between two
** specific times.  The buffer is of type ADSMarkers, pointing at an array
** of markers.  The caller requests a number of markers, the routine gives
** back the number it actually transferred, or -ve error code.
**
** fh           The handle for the file
** chan         The data channel to get times from
** markP        Where to store the times we find
** max          The maximum number of events to collect
** sTime,eTime  The start and end of the period we want data in
** pFiltMask     Marker filter, or NULL.
*/
SONAPI(long) SONGetMarkData(short fh, WORD chan, TpMarker markP, long max,
        TSTime sTime, TSTime eTime, TpFilterMask pFiltMask)
{
    long block;                 /* used to hold next block */
    long dataLeft;              /* data poins left in a block */
    long count = 0;             /* number of points returned */
    long nBytes;                /* size of a data item */
    TpMarker markRdP;           /* marker read pointer */
    TpMarker markWrP = markP;   /* writing pointer WHY COPY IT? */
    TpDataBlock bufP;           /* local copy of buffer pointer */
    TDataKind kind;             /* will hold the channel type */
    TpChannel cP;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(max > 0);
    
    bufP = files[fh].bufferP;                   /* get local pointers */
    cP = SONChanPnt(fh, chan);
    kind = cP->kind;                          /* get the channel type */

    if ((kind != Marker) &&
        (kind != AdcMark) &&
        (kind != RealMark) &&
        (kind != TextMark))
        return SON_NO_CHANNEL;      /* was not a marker */

    block = SONFindBlock(fh, chan, sTime, eTime);  /* find a block */
    if (block <= 0)
        return block;                   /* return with error if none */

    nBytes = SONItemSize(fh,chan);

    do                  /* loop once per data buffer */
    {
        long iErr = SONReadBlock(fh,chan, block); /* read block to bufferP */
        if (iErr < 0)
            return iErr;                    /* error if we fail */
        dataLeft = bufP->items;             /* number of items in the buffer */
        markRdP = bufP->data.markData;      /* set data pointer */

        if ((count == 0) && dataLeft)       /* if first block */
            while (markRdP->mark < sTime)   /* while time not in range...*/
            {
                markRdP = (TpMarker)movePtr(markRdP,nBytes); /* ...move pointer &... */
                dataLeft--;                 /* ...reduce data count */
            }

        while (dataLeft && (count<max))     /* loop to collect data */
        {
            if (markRdP->mark > eTime)      /* if beyond wanted time...*/
                return count;               /* ...return count of events */
            if ( !pFiltMask ||              /* is filter defined */
                  SONFilter(markRdP, pFiltMask))   /* do we want this data */
            {
                *markWrP++ = *markRdP;      /* add to output */
                count++;
            }
            markRdP=(TpMarker)movePtr(markRdP,nBytes);
            dataLeft--;
        }

        block = bufP->succBlock;            /* next block to read */
        if (block > cP->lastBlock)          /* for GRB family */
            block = CHAINEND;               /* as chain points at vacuum! */
    }
    while( (count < max)                    /* buffer not full and...*/
           && (block > 0)                   /* another block available and*/
           && (bufP->endTime <= eTime) );   /* not past the end time */
    return count;
}

/****************** S O N G e t A D C D a t a ( ) ********************
**
** SONGetADCData:- Used to collect ADC data into a buffer between two
** specific times.  The buffer is of type ADSAdc, pointing at an array
** of data.  The caller requests a number of samples, the routine gives
** back the number it actually transferred, or -ve error code.
**
** fh           The handle for the file
** chan         The data channel to get times from
** adcDataP^    Where to store the times we find
** max          The maximum number of events to collect
** sTime,eTime  The start and end of the period we want data in
** pbTime       Point to val set to actual time of first ADC point in buff
** pFiltMask     Pointer to marker filter, or NULL.
*/
SONAPI(long) SONGetADCData(short fh, WORD chan, TpAdc adcDataP, long max,
                  TSTime sTime, TSTime eTime, TpSTime pbTime, 
                  TpFilterMask pFiltMask)
{
    long block;         /* to find the data in the file */
    TSTime adcTime;     /* time per adc sample */
    TSTime nextTime;    /* predicted time of next data block */
    long dataIndex;     /* index into the data buffer to the required point */
    long lastIndex;     /* index of last item in this block we want */
    long words;         /* number of items to copy */
    long count = 0;     /* number of items copied */
    long spaceIndex;    /* max index into array for which there is space */
    long lErr;
    TpChannel cP;       /* Our chan and datablock pointers */
    TpDataBlock bufP;
    TDataKind kind;     /* type of the channel */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(max > 0);
    assert(pbTime != NULL);
    
    kind = SONChanKind(fh,chan);
    if ((kind != Adc) && (kind != AdcMark))
        return SON_NO_CHANNEL;
    if ((block = SONFindBlock(fh, chan, sTime, eTime))<=0)
        return block;
    lErr = SONReadBlock(fh, chan, block);
    if ( lErr < 0 )
        return lErr;
    cP = SONChanPnt(fh, chan);
    bufP = files[fh].bufferP;           /* pointer to the data buffer */

    adcTime = SONChanDivide(fh, chan);  /* time per Adc sample */

    /*
    ** Code to deal with ADC marker data. We must allow for overlapped markers
    ** and return contiguous data that spans several markers if we can.
    */
    if (SONChanKind(fh,chan) == AdcMark)
    {
        int nBytes = SONItemSize(fh,chan);      /* Number of bytes per item */
        int nPoints = (WORD)(cP->nExtra >> 1);  /* Number of ADC points per item */
        TSTime adcMkSTime = sTime - ((nPoints-1)*adcTime);/* Start time for search */
#ifdef USEHANDLES
        THandle hAdcMarkP;  /* and handle as necessary */
        TpAdcMark adcMarkP = (TpAdcMark ) HAlloc( nBytes, &hAdcMarkP );
#else
        TpAdcMark adcMarkP = (TpAdcMark ) F_malloc( nBytes );
#endif
        if ( adcMarkP == NULL )
            return (SON_OUT_OF_MEMORY);

        while ((adcMkSTime < eTime) &&      /* might still get data, and... */
               (count < max))               /* ...buffer not full */
        {
            TSTime lTime;                   /* time of this WaveMark */
            long lStOff = 0;                /* index of first useful point */
            long lLastIndex;                /* last useable index in the buffer */
            long lCopy;                     /* number of points to add */
            TSTime lDataTime;               /* first useable point time */

            /* get the next marker that is possible */
            lErr = SONGetExtMarkData(fh, chan, (TpMarker)adcMarkP, 1, adcMkSTime,
                                                            eTime, pFiltMask);
            if (lErr <= 0)                  /* no more data found */
            {
                if (lErr < 0)               /* negative return is an error */
                    count = lErr;           /* if error, return it */
                break;                      /* get out of the while loop */
            }

            lTime = adcMarkP->m.mark;       /* Time of marker */
            if (lTime < sTime)              /* first point possibly not wanted */
                lStOff = ((sTime - lTime + adcTime -1) / adcTime);

            lDataTime = lTime + adcTime * lStOff;  /* first data point time */
            if (count == 0)                 /* is this the first data point? */
                *pbTime = lDataTime;        /* set as first time in buffer */
            else
            {
                if (lDataTime != adcMkSTime)
                    break;                  /* not contiguous, so give up */
            }

            /* There may be data beyond last wanted point */
            lLastIndex = (eTime - lTime) / adcTime;     /* last point we want */
            if (lLastIndex >= nPoints)      /* limit the index to actual data */
                lLastIndex = nPoints-1;     /* use long, eTime could be large */

            /* Now work out how many points we could add to buffer */
            lCopy = lLastIndex - lStOff + 1;
            if (lCopy + count > max)        /* make sure not too many */
                lCopy = max - count;

            if (lCopy > 0)                  /* anything to copy? */
            {
                F_memcpy(adcDataP + count, adcMarkP->a + lStOff, lCopy * 2);
                count += lCopy;
                adcMkSTime = lDataTime + lCopy * adcTime; /* next point */
            }
            else
                break;                      /* must give up if no points */
        }
#ifdef USEHANDLES
        HFree(hAdcMarkP);   /* Free allocated memory */
#else
        F_free(adcMarkP);
#endif
        return count;       /* Return with ADC data */
    }

/*
** Now back to normal Adc data; calculate first useful value...
** there must be one as SONFindBlock has data.
*/
    if (bufP->startTime>sTime)
        dataIndex = 0;      /* start with the first point */
    else
        dataIndex = ((sTime-bufP->startTime+adcTime-1)/adcTime);

    *pbTime = bufP->startTime + (adcTime*dataIndex); /*first point time*/

    do
    {
        if (bufP->endTime >= eTime)     /* do we finish in this block? */
            lastIndex = ((eTime-bufP->startTime)/adcTime);
        else
            lastIndex = (bufP->items - 1);
        /*
        ** Now check we have space for these items.
        */
        spaceIndex = (max - count + dataIndex - 1);  /* max index */
        if (lastIndex > spaceIndex)
            lastIndex = spaceIndex;
        words = (lastIndex + 1 - dataIndex);  /* number of words to copy */
        F_memcpy(adcDataP+count, bufP->data.int2Data+dataIndex, words*2);
        count += words;                     /* increase count done */
        block = bufP->succBlock;            /* the next block */
        if (block > cP->lastBlock)          /* for GRB family */
            block = CHAINEND;

        if ((bufP->endTime>=eTime)          /* if next buffer too late */
            || (count == max)               /* or buffer is full */
            || (block<0))                   /* or no next block */
            return count;                   /* job done */

        nextTime = bufP->endTime + adcTime; /* next block if contiguous */
        lErr = SONReadBlock(fh, chan, block);  /* read block */
        if (lErr < 0)                       /* report any errors */
            return lErr;
                               
        dataIndex = 0;                      /* next index is buffer start */
    }while( bufP->startTime == nextTime );  /* while contiguous */

    return count;
}

/******************* S O N G e t E x t M a r k D a t a ( ) ****************
**
** SONGetExtMarkData:- Used to collect ADC, Real or TextMarker data into
** a buffer. Only data between two specific times is returned.  The buffer
** is of type *TMarker, pointing at an array of data.  The caller requests
** a number of samples, the routine gives back the number it actually
** transferred, or -ve error code.
**
** fh           The handle for the file
** chan         The data channel to get data from
** markDataP    A pointer to where to store the data we find
** max          The maximum number of items to collect
** sTime,eTime  The start and end of the period we want data in
** pFiltMask     Pointer to marker filter or NULL.
*/
SONAPI(long) SONGetExtMarkData(short fh, WORD chan, TpMarker markDataP,
                      long max, TSTime sTime, TSTime eTime, 
                      TpFilterMask pFiltMask)
{
    long block;         /* to find the data in the file */
    long count = 0;     /* number of items copied */
    TDataKind kind;     /* the channel type */
    WORD nBytes;        /* size of the item in bytes */
    TpDataBlock bufP;   /* points to data buffer */
    TpMarker markRdP;   /* data pointer for reading from buffer*/
    TpChannel cP;

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    assert(max > 0);
    
    bufP = files[fh].bufferP;                      /* get local pointers */
    cP = SONChanPnt(fh, chan);
    kind = cP->kind;
    if ((kind != AdcMark) &&
        (kind != RealMark) &&
        (kind != TextMark) &&
        (kind != Marker))
        return SON_NO_CHANNEL;

    if ( (block = SONFindBlock(fh, chan, sTime, eTime)) <= 0 )
        return (short) block;

    nBytes = SONItemSize(fh, chan);     /* Number of bytes in item */

    do
    {   /* loop once per data buffer */
        long dataLeft;
        short err;

        if ((err = SONReadBlock(fh, chan, block))<0) /* fill data buffer */
            return err;

        markRdP = bufP->data.markData;  /* set read pointer position */
        dataLeft = bufP->items;         /* number of items in the buffer */

        if ( (count==0) && dataLeft )   /* first time round, skip to data */
            while (markRdP->mark < sTime)   /* find first wanted item */
            {
                markRdP = (TpMarker)movePtr(markRdP,nBytes);
                dataLeft--;
            }

        while (dataLeft && (count<max))     /* loop to collect data */
        {
            if (markRdP->mark > eTime)      /* are we past wanted data? */
                return count;               /* if so, return events */

            if ( !pFiltMask ||               /* no filter defined? */
                 SONFilter(markRdP, pFiltMask))    /* or defined and true */
            {
                F_memcpy(markDataP, markRdP, nBytes);  /* add to output */
                markDataP=(TpMarker)movePtr(markDataP,nBytes); /* on to next marker */
                count++;                    /* increase count of data */
            }
            dataLeft--;
            markRdP=(TpMarker)movePtr(markRdP,nBytes);/* on to the next marker */
        }

        block = bufP->succBlock;            /* next data block */
        if (block > cP->lastBlock)    /* for GRB family */
            block = CHAINEND;
    }
    while( (count < max)                    /* buffer not yet full */
           && (block > 0)                   /* and more on disk */
           && (bufP->endTime <= eTime) );   /* and not out of time */

    return count;
}/* GetExtMarkData */


/******************* S O N S e t M a r k e r ( ) ******************
**
** SONSetMarker:- Used to replace an existing marker with new data. The
** marker at the specified time is found and it's data replaced with
** the new data supplied. If the time for the new marker is different to
** the old, the disk block is sorted to keep markers in time order. The
** new time may be altered slightly to prevent two markers at one time.
** If the new time does not fit within the current block, an error is
** returned. Returns 1 if marker replaced, zero if marker with correct
** time was not found or new time was too different, or a -ve error code.
**
** fh           The handle for the file
** chan         The data channel to get data from
** time         The time for the marker to be changed
** markP        A pointer to the replacement marker data
*/
SONAPI(short) SONSetMarker(short fh, WORD chan, TSTime time, 
                                                TpMarker markP, WORD size)
{
    long        block;          /* offset to the data in the file */
    TDataKind   kind;           /* the channel type */
    long        nBytes;         /* size of the item in bytes */
    TpDataBlock bufP;           /* points to data buffer */
    TpMarker    markRdP;        /* data pointer for reading from buffer*/
    TpChannel   cP;             /* Channel info pointer */
    short       err;            /* GP error code */
    WORD        item;           /* GP counter */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    bufP = files[fh].bufferP;                      /* get local pointers */
    cP = SONChanPnt(fh, chan);
    kind = cP->kind;
    if ((kind != AdcMark) &&
        (kind != RealMark) &&
        (kind != TextMark) &&
        (kind != Marker))
        return SON_NO_CHANNEL;
    nBytes = SONItemSize(fh, chan);              /* Number of bytes in item */
    if (size > (WORD)nBytes)        /* Must not be replacing more than this */
        return SON_NO_CHANNEL;
    if (size < 4)
        return SON_NO_CHANNEL;         /* Must have at least 4 bytes to replace */

    if ( (block = SONFindBlock(fh, chan, time, time)) <= 0 )
        return (short) block;         /* Find block covering specified time */
    err = SONReadBlock(fh, chan, block);      /* fill data buffer */
    if (err < 0)
        return err;

    item = 1;               /* Find the item in the block with correct time */
    markRdP = bufP->data.markData;            /* Starting at the first item */
    while ((markRdP->mark < time) && (item < bufP->items))
    {
        item++;                                     /* move on to next item */
        markRdP = (TpMarker)movePtr(markRdP, nBytes);
    }
    if ((item > bufP->items) || (markRdP->mark != time))
        return 0;                  /* Return if found nothing or wrong time */

    if (item == 1)               /* Check new time doesn't overlap previous */
    {                                /* Look at previous block if necessary */
        if (bufP->predBlock != CHAINEND)    /* but not if at start of chain */
        {
            err = SONReadBlock(fh, chan, bufP->predBlock);
            if (err < 0)
                return err;
            if (bufP->endTime >= markP->mark)        /* must be after predBlock */
                return 0;
            if ((err = SONReadBlock(fh, chan, block)) < 0)
                return err;
        }
    }
    else
    {
        TpMarker temp = (TpMarker)movePtr(markRdP, -nBytes);  /* Point to prev marker */

        if (markP->mark <= temp->mark)               /* and check it's time */
            return 0;
    }

    if (item == bufP->items)     /* and also check it is OK for next marker */
    {                                    /* Look at next block if necessary */
        if ((bufP->succBlock != CHAINEND) &&   /* but not if at start of chain */
                (bufP->succBlock <= cP->lastBlock))       /* for GRB family */
        {
            err = SONReadBlock(fh, chan, bufP->succBlock);
            if (err < 0)
                return err;
            if (bufP->endTime <= markP->mark)        /* must be after predBlock */
                return 0;
            err = SONReadBlock(fh, chan, block);
            if (err < 0)
                return err;
        }
    }
    else
    {
        TpMarker temp = (TpMarker)movePtr(markRdP, nBytes);   /* Point to next marker */

        if (markP->mark >= temp->mark)                   /* check it's time */
            return 0;
    }

    F_memcpy(markRdP, markP, size);          /* copy new data into position */
    err = SONWrite(fh, bufP, cP->phySz, block);   /* Write back to the file */
    if (err < 0)
        return err;

    return 1;                             /* Return indicating job was done */
}

/***************** S O N G e t V e r s i o n ( ) ********************
**
** return the file version number or -1 if no file
*/
SONAPI(int) SONGetVersion(short fh)
{
    assert((unsigned)fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    if (!files[fh].opened)          /* check we have a file */
        return -1;                  /* no file, so -1 */
    else
        return (int)files[fh].headP->systemID;
}

/***************** S O N G e t E x t r a D a t a S i z e ( ) ********
**
** return the number of bytes of extra data in the file
*/
SONAPI(long) SONGetExtraDataSize(short fh)
{
    assert((unsigned)fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);
    if (!files[fh].opened)          /* check we have a file */
        return 0;                  /* no file, so 0 */
    else
        return (long)files[fh].headP->extraData;
}

/***************** S O N G e t E x t r a D a t a ( ) ****************
**
** This transfers the extra data (if any) between the data file and the
** nominated buffer.
*/
SONAPI(short) SONGetExtraData(short fh, void FAR *buff, WORD bytes,
                                               WORD offset, BOOLEAN writeIt)
{
    TpFileHead fileHeadP;                        /* local head pointer */
    long position;                          /* space for file position */

    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (!files[fh].opened)                  /* check we have a file */
        return SON_NO_FILE;                     /* must have a file */
    fileHeadP = files[fh].headP;
    if (!fileHeadP->extraData               /* no extra data or... */
        || ((WORD)(bytes+offset)>fileHeadP->extraData)) /* ...out of range */
        return SON_NO_EXTRA;

    position = (long)offset + DISKBLOCK + CHANSIZE; /* file position */
    if (writeIt)                            /* see if read or write */
        return SONWrite(fh, buff, bytes, position);
    else
        return SONRead(fh, buff, bytes, position);
}

/**************** S O N C h a n D e l e t e ( ) ************************
**
** mark the selected channel as deleted if it was in use and had storage
** space associated with it.  This allows us to junk a channel without
** releasing the space to DOS.  The channel is then free for use (but it
** must be used with the same size disk blocks please.
*/
SONAPI(short) SONChanDelete(short fh, WORD chan)
{
    TpChannel   cP;                 /* Pointer to channel data */
    short       err;                /* GP error flag */

    assert(chan<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    if (chan >= SONMAXCHANS)    /* ignore silly request */
        return 0;
    cP = SONChanPnt(fh, chan);  /* get a pointer to the channel */
    if (cP->kind == ChanOff)
        return 0;               /* nothing to do if off already */
    cP->kind = ChanOff;         /* say channel is off */
    if (cP->blocks > 0)         /* if there are existing blocks */
    {
        if (cP->delSize != 0)   /* if previous deleted chain */
        {             /* Append old deleted chain to end of current data */
            err = SONSetSucc(fh, cP->lastBlock, cP->nextDelBlock);
            if (err < 0)        /* If we failed, give up right now ! */
                return err;
        }
        cP->delSize += cP->blocks;   /* update size of deleted chain */
        cP->nextDelBlock = cP->firstBlock;     /* set start of chain */
    }
    cP->blocks = 0;                        /* no blocks left */
    cP->firstBlock = CHAINEND;
    cP->lastBlock = CHAINEND;
    files[fh].speedP[chan].prevBlock = CHAINEND;
    SONUpdateMaxTimes(fh);
    return SONUpdateStart(fh);  /* write changes to the file */
}


/************************ S O N L a s t T i m e ***************************
**
** Find the last item on the channel that occured before the set time.
**
** fh       The son file handle to use
** channel  The data channel to search
** sTime    The time to start the search at
** eTime    Search upto and including this time (NOTE: eTime < sTime !)
** psVal    Pointer to the event level for EventBoth data, pointer to
**          adc data value for Adc channel. Undefined for other channels.
** mbyte    The marker bytes
** pbMk     Returned TRUE if it is a marker channel
** pFiltMask Pointer to marker filter, or NULL to accept all.
**
** Returns the time if a data item is found, or -1 if not found or an error
*/

SONAPI(TSTime) SONLastTime(short fh, WORD channel, TSTime sTime, TSTime eTime,
                    short FAR *psVal, TMarkBytes FAR *mBytes, 
                    TpBOOL pbMk, TpFilterMask pFiltMask)
{
    long block;             /* The disk position to read for the data */
    TDataKind kind;         /* what kind of channel doya call this? */
    long theTime = -1;      /* assume that we won't find anything */
    BOOLEAN isMarker = FALSE;

    assert(channel<SONMAXCHANS);
    assert(fh >= 0);
    assert(fh < MAXFILES);
    assert(files[fh].defined);
    assert(files[fh].opened);

    kind = SONChanKind(fh, channel);/*get channel type*/
    if (kind == ChanOff)  /* FindBlock tells us where to read from */
        return -1;

    *pbMk = isMarker = (BOOLEAN)((kind == AdcMark) || (kind == Marker) ||
                                (kind == RealMark) || (kind == TextMark));

    /*
    ** first we must find which block the data is in so we can ask for
    **  the correct time from the channel
    */
    
    block = SONFindBlock(fh, channel, sTime, LONG_MAX);  /* find the time */
    if (block <= 0)         /* fail means no data, or past all data */
        block = SONChanPnt(fh, channel)->lastBlock; /* So try last block */

    if (block > 0)                      /* found a useful block*/
    {
        TpDataBlock pDB = files[fh].bufferP;    /* local pointer to data */
        if (SONReadBlock(fh, channel, block) != 0)  /* read it */
            return -1;                  /* give up if could not read */

        if (pDB->startTime >= sTime)    /* is our data in this block?*/
        {
            block = pDB->predBlock;     /* NO! we want the previous one*/

            if (block <= 0)             /* if no block, then return */
                return -1;

            if (SONReadBlock(fh, channel, block) != 0)  /* this must have it*/
                return -1;
        }

        /*
        ** Now, the data we require is in the block we just read unless
        ** this is a marker and we have a filter set (in which case we
        ** must search backwards for the data we want).
        */
        if (pDB->startTime < sTime)     /* check must be a value */
        {
            if (kind == Adc)            /* just work out time for Adc */
            {
                TSTime adcTime = SONChanDivide(fh, channel);    /* interval */
                int adcIndex = (int)((sTime - 1 - pDB->startTime) / adcTime);
                if (adcIndex<0)         /* make sure no underflow */
                    adcIndex = 0;
                else
                    if ((WORD)adcIndex>=pDB->items)
                        adcIndex=(int)pDB->items-1;
                theTime = pDB->startTime + adcIndex * adcTime;
                *psVal = pDB->data.int2Data[adcIndex]; /*starts at 1*/
            }
            else
            {
                WORD skipBytes = SONItemSize(fh, channel);
                WORD count = pDB->items;
                TpMarker srcPtr = (TpMarker)((LPSTR)(&pDB->data) +
                                                   (count - 1) * skipBytes);

                while (srcPtr->mark >= sTime)   /*found it?*/
                {
                    srcPtr = (TpMarker)((LPSTR)srcPtr - skipBytes);
                    count--;
                }

                if (isMarker && (pFiltMask != NULL))
                {
                    while ((srcPtr->mark >= eTime) &&
                                    !SONFilter(srcPtr, pFiltMask))
                    {
                        srcPtr = (TpMarker)((LPSTR)srcPtr - skipBytes);
                        if (--count == 0)   /* reached buffer start? */
                        {
                            block = pDB->predBlock; /* try previous block */
                            if (block > 0)
                            {
                                if (SONReadBlock(fh, channel, block) != 0)
                                    return -1; /* give up if could not read */
                                count = pDB->items;
                                srcPtr = (TpMarker)((LPSTR)(&pDB->data) +
                                                   (count - 1) * skipBytes);
                            }
                            else
                                return -1;  /* nothing to find */
                        }
                    }
                }

                theTime = srcPtr->mark;

                if (kind == EventBoth)
                    *psVal = (short)(((pDB->chanNumber & 256) != 0) ^ 
                                                    (count & 1));

                if (isMarker)
                {
                    (*mBytes)[0] = srcPtr->mvals[0];
                    (*mBytes)[1] = srcPtr->mvals[1];
                    (*mBytes)[2] = srcPtr->mvals[2];
                    (*mBytes)[3] = srcPtr->mvals[3];
                }
            }              /* if earlier than search range return not found */
            return ((theTime >= eTime) ? theTime : -1);
        }
    }
    return -1;
}

/* 
**
** Code to handle filters for markers for SON filing system.
**
** You can define TFilterElt to be unsigned byte, word or long if there
** is any computational advantage on the target computer. If you set
** word, ELMSK = 0xf, ELSHF=4. If long, ELMSK=0x1f, ELSHF=5.
*/

#define MAXITEM 255                 /* maximum item in a mask */
#define ELMSK 0x7                   /* mask for the element */
#define ELSHF 3                     /* shift to find element number */
#define MASKSZ 32                   /* number of TFilterElt in mask */

/*********************************************************************
**
** S O N f i l t e r ( )
**
** Used to decide if a given marker is part of a set of markers. This
** uses the cAllSet flag to test if the layer needs testing. As layers
** 1-3 will in general not need checking, this should save quite a bit
** of calculation. We drop out as soon as a layer fails the test.
**
*/
SONAPI(int) SONFilter(TpMarker pM, TpFilterMask pFM)
{
    int i;              /* used to keep track of each layer */
    for (i=0; i<4; i++) /* try each layer in turn */
        if (pFM->cAllSet[i] == 0)         /* only test if set not full */
        {
            int index = ((unsigned char)pM->mvals[i]) >> ELSHF; /* element index in table */
            TFilterElt mask = (TFilterElt)(1 << (pM->mvals[i] & ELMSK));
            if ((mask & pFM->aMask[i][index]) == 0) /* not in? */
                return 0;                   /* if not in, give up */
        }                                   /* if not found, give up */
    return 1;           /* if we get here, the marker is included */
}

/*********************************************************************
**
** t e s t A l l S e t ( )
**
** This function checks if an entire layer is set. If so, it sets the
** layer flag to allow the filter function to save time. layer is taken
** as being already tested.
*/
static void testAllSet(TpFilterMask pFM, int layer)
{
    int i;
    pFM->cAllSet[layer] = 0;        /* assume the worst */
    for (i=0; i<MASKSZ; i++)        /* test for full bitmap */
        if (pFM->aMask[layer][i] != (TFilterElt)0xff)
            return;                 /* if any unset, we quit */
    pFM->cAllSet[layer] = 1;        /* all are set if we get here */
}

/*********************************************************************
**
** d o L a y e r ( )
**
** This function operates on an entire layer to set, clear or invert it.
**
** pFM   Pointer to the mask bitmap.
** layer The layer to operate on (0-3). It is assumed that the layer does
**       not exceed this range.
** set   This defines the operation: SON_FCLEAR clears the layer,
**       SON_FSET sets the layer and SON_FINVERT inverts it.
**       -ve value reads the layer, returning 1 if whole layer set, else 0.
**
** As this function would not be used in any kind of a speed critical way,
** it is coded to be efficient in space rather than efficient in time.
*/
static int doLayer(TpFilterMask pFM, int layer, int set)
{
    if (set < 0)
        return pFM->cAllSet[layer];
    else
    {
        int i;                                  /* loop counter */
        TFilterElt FAR *pE = &pFM->aMask[layer][0]; /* point at first element */
        for (i=0; i<MASKSZ; i++, pE++)          /* run round all elements */
            switch (set)
            {
            case SON_FCLEAR:
                *pE = (TFilterElt)0;
                break;
    
            case SON_FSET:
                *pE = (TFilterElt)0xff;
                break;
    
            case SON_FINVERT:
                *pE ^= (TFilterElt)0xff;
                break;
            }
        testAllSet(pFM, layer);                 /* see if all set */
        return 0;
    }
}

/*********************************************************************
**
** d o I t e m ( )
**
** This function operates on one item of a layer to set, clear or invert it.
**
** pFM   Pointer to the mask bitmap.
** layer The layer to operate on (0-3). It is assumed that the layer does
**       not exceed this range.
** item  The item in the layer (0-255). It is assumed that the item does
**       not exceed this range.
** set   This defines the operation: SON_FCLEAR clears the layer,
**       SON_FSET sets the layer and SON_FINVERT inverts it.
**       -ve value reads the item, returning 1 if set, else 0.
**
** As this function would not be used in any kind of a speed critical way,
** it is coded to be efficient in space rather than efficient in time.
*/
static int doItem(TpFilterMask pFM, int layer, int item, int set)
{
    int index = item >> ELSHF;  /* calculate index and bitmask */
    TFilterElt mask = (TFilterElt)(1 << (item & ELMSK));
    TFilterElt FAR *pE  = & pFM->aMask[layer][index];    /* element to use */

    if (set < 0)            /* is this a read operation? */
        return ((*pE & mask) != (TFilterElt)0);
    else
    {
        switch (set)            /* must be clear, set or invert */
        {
        case SON_FCLEAR:        /* clear the bit */
            *pE &= ~mask;
            break;
    
        case SON_FSET:          /* set the bit */
            *pE |= mask;
            break;
    
        case SON_FINVERT:       /* invert the bit */
            *pE ^= mask;
            break;
        }
        testAllSet(pFM, layer); /* set layer flag */
        return 0;
    }
}

/*******************************************************************
**
** SONFControl will set/clear individual bits in one or all layers
**
** pFM      Pointer to the structure to use
** layer    The layer to select: 0-3. -ve means all layers
** item     The item in the layer, 0-255. -ve means all, >255 is an error.
** set      -ve = read the item, 0=clear, 1=set, 2=invert. >2 is an error.
**
** If set is -ve, we read data for the layer (0-3) and the item (0-255)
** and return 0 or 1. If item is out of range we return -1. If ask for whole
** layer or all layers then return 1 if all items are set else 0.
** 
** If set is 0 or >0, we reset or set the bit (or layer if item is -ve)
** and return 0. If there is a problem, return -1.
**
** As speed is not critical here, we write the code to be space efficient.
*/
SONAPI(int) SONFControl(TpFilterMask pFM, int layer, int item, int set)
{
    if ((layer > 3) ||          /* give up if stupid layer found, or... */
        (item > MAXITEM) ||     /* ...a stupid item in the layer, or... */
        (set > SON_FINVERT))    /* ...an unknown operation */
        return -1;              /* and return an error */

    if (layer < 0)              /* is this all layers? */
    {
        int i;
        int ret = 1;
        for (i=0; i<4; i++)     /* for all layers */
        {
            if (item < 0)       /* whole layer */
                ret &= doLayer(pFM, i, set);
            else
                ret &= doItem(pFM, i, item, set);
        }
        return ret;
    }
    else 
    {
        if (item < 0)          /* whole layer operation? */
            return doLayer(pFM, layer, set);
        else
            return doItem(pFM, layer, item, set);
    }
    
    return 0;
}

/*******************************************************************
**
** SONFEqual returns TRUE if the 2 filter masks are equal, else FALSE
**
*/
SONAPI(BOOLEAN) SONFEqual(TpFilterMask pFiltMask1, TpFilterMask pFiltMask2)
{
    short layer;
    short elt;
    
    for (layer = 0; layer < 4; layer++)
    {
        if (pFiltMask1->cAllSet[layer] != pFiltMask2->cAllSet[layer])
            return FALSE;
    }

    for (layer = 0; layer < 4; layer++)
    {
        for (elt = 0; elt < 32; elt++)
        {
            if (pFiltMask1->aMask[layer][elt] != 
                            pFiltMask2->aMask[layer][elt])
                return FALSE;
        }
    }
    
    return TRUE;    /* if get here then filter masks are equal */
}
