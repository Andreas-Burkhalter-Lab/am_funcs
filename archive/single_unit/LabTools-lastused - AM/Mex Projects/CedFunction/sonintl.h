/*
** sonintl.h
**
******************************************************************************
**
**  Definitions of the file structures and other internals for the SON
**  filing system. This is not normally required by library users, but
**  is supplied for specialised use. The file machine.h (included in son.h)
**  provides some common definitions across separate platforms.
**
**  This should be included AFTER son.h.
**
*/

#ifndef __SON__
#include "son.h"
#endif

#ifndef __SONINTL__
#define __SONINTL__

#define LSTRING(size) union{unsigned char len;char string[size+1];}
#define REVISION 5
#define MAXFILES 32                  /* The maximum number of files allowed */
#define DISKBLOCK 512
#define ROUND_TO_DB(num) (((num)+DISKBLOCK-1)&0xfe00)
#define LENCOPYRIGHT  10
#define LENSERIALNUM   8
#define COPYRIGHT "(C) CED 87"

/* the field osFormat is a 16 bit word, used to detect file format, so byte */
/* swapping will not affect this value */
#define DOSFORMAT   0x0000
#define MACFORMAT   0x0101
#if defined(macintosh) || defined(_MAC)
#define OSFORMAT MACFORMAT
#else
#define OSFORMAT DOSFORMAT
#endif

#define CHAINEND -1  /* The value that indicates the end of chain of blocks */

typedef LSTRING(SON_CHANCOMSZ) TChanComm;
typedef LSTRING(SON_COMMENTSZ) TComment;
typedef LSTRING(SON_TITLESZ) TTitle;
typedef TComment TFileComment[SON_NUMFILECOMMENTS];  /* file comment for header */
typedef LSTRING(SON_UNITSZ) TUnits;            /* units string for adc channels */

/*
** Macro to mode a generic pointer on by n bytes
*/
#define movePtr(p,n) ((void FAR *)((TpStr)(p) + (n)))

#if defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)
#pragma pack(2)
#endif

/*
** Now a structure which defines the first disk block of a file.  We
** pad the structure out to 512 bytes as we assume that all files are
** efficient if we make reads and writes multiples of 512
** 2.5.91 extra field osFormat added, and padding reduced to 72 chars
*/
typedef struct      /* first disk block of file */
{
    short systemID;               /* filing system revision level */
    char copyright[LENCOPYRIGHT]; /* space for "(C) CED 87" */
    char serialNum[LENSERIALNUM]; /* space for serial number */
    WORD usPerTime;               /* microsecs per time unit */
    WORD timePerADC;              /* time units per ADC interrupt */
    short fileState;              /* condition of the file */
    long firstData;               /* offset to first data block */
    short channels;               /* maximum number of channels */
    WORD chanSize;                /* memory size to hold chans */
    WORD extraData;               /* No of bytes of extra data in file */
    WORD bufferSz;                /* Not used on disk; bufferP in bytes */
    WORD osFormat ;               /* either 0x0101 for Mac, or 0x00 for PC */
    TSTime maxFTime;              /* max time in the data file */
    char pad[68];                 /* padding for the future */
    TFileComment fileComment;     /* what user thinks of it so far */
} TFileHead;

typedef TFileHead FAR * TpFileHead;

/*
** TChannel is a structure which tells us about an individual channel
** of data.  An array of these follows the header block of the file.
** Note that the look up table only exists in memory.  It is meaning
** less on disk, unless we actually store it on disk at some future
** time.
*/
typedef struct
{
    WORD delSize;       /* number of blocks in deleted chain, 0=none */
    long nextDelBlock;  /* if deleted, first block in chain pointer */
    long firstBlock;    /* points at first block in file */
    long lastBlock;     /* points at last block in file */
    WORD blocks;        /* number of blocks in file holding data */
    WORD nExtra;        /* Number of extra bytes attached to marker */
    short preTrig;      /* Pre-trig points for ADC Marker data */
    short free0;        /* Keeps space OK */
    WORD phySz;         /* physical size of block written =n*512 */
    WORD maxData;       /* maximum number of data items in block */
    TChanComm comment;  /* string commenting on this data */
    long maxChanTime;   /* last time on this channel */
    long free1;         /* spare space for the future, set to 0 */
    short phyChan;      /* physical channel used */
    TTitle title;       /* user name for channel */
    float idealRate;    /* ideal rate:ADC, estimate:event */
    TDataKind kind;     /* data type in the channel */
    unsigned char pad;  /* padding just to keep up with Pascal... */

    union               /* now section which changes with the data */
    {
        struct
        {                       /* Data for ADC and ADCMark channels */
            float scale;
            float offset;       /* to convert to units */
            TUnits units;       /* channel units */
            WORD divide;        /* from ADC int rate */
        } adc;
        struct
        {                       /* only used by EventBoth channels */
            BOOLEAN initLow;    /* initial event state */
            BOOLEAN nextLow;    /* expected state of next write */
        } event;
        struct
        {                       /* This one for real marker data */
            float min;          /* expected minimum value */
            float max;          /* expected maximum value */
            TUnits units;       /* channel units */
        } real;                 /* NB this is laid out as for adc data */
    } v;

} TChannel;

typedef TChannel FAR * TpChannel;

/*
** Now a structure, being all the channels in an array.  This is saved
** on disk starting at offset 512 into the file. We also define a Macro for
** the size of this structure, rounded up to a disk block. Apologies for
** the 0xfe00, but this is simplest way to express the result.
*/
typedef TChannel TChannels[SONMAXCHANS];
#define CHANSIZE ((sizeof(TChannels) + DISKBLOCK -1) & 0xFE00)

/*
** The data is stored in blocks (again multiples of 512 bytes long)
** on disk.  All the blocks have an identical header, but the rest
** depends on what the data is.
**
** On the Mac you're not allowed structs > 32k, so the original definition of
** TDataBlock will not compile. So I have defined ADCdataBlkSize etc to be
** half the size
*/
#if defined(macintosh) || defined(_MAC)
#define ADCdataBlkSize  16000
#define timeDataBlkSize 8000
#define markDataBlkSize 4000
#else
#define ADCdataBlkSize  32000
#define timeDataBlkSize 16000
#define markDataBlkSize 8000
#endif


typedef struct
{
    long   predBlock;     /* predecessor block in the file */
    long   succBlock;     /* following block in the file */
    TSTime startTime;     /* first time in the file */
    TSTime endTime;       /* last time in the block */
    WORD   chanNumber;    /* The channel number in the block */
    WORD   items;         /* Actual number of data items found */
    union
    {
        TAdc      int2Data [ADCdataBlkSize];    /* ADC data */
        TSTime    int4Data [timeDataBlkSize];   /* time data */
        TMarker   markData [markDataBlkSize];   /* marker data */
        TAdcMark  adcMarkData;                  /* ADC marker data */
    } data ;
} TDataBlock;

typedef TDataBlock FAR * TpDataBlock;
#define SONDBHEADSZ 20  /* size of the header data for each block */

typedef struct      /* this structure is used to speed up disk reading */
{
    long prevBlock;
    long prevPred;
    long prevSucc;
    long prevStart;
    long prevEnd;
} TSpeedPtr;

typedef TSpeedPtr FAR * TpSpeedPtr;

/*
** Structure used to define a file for us to actually use
** we use a DOS file handle to refer to the file.
*/

#ifdef USEHANDLES        /* Types required for memory handle work */
                         /* Used to define THandle here - now machine.h */
#endif                   /* End of handle type definitions */

typedef struct
{
    BOOLEAN opened;     /* set true if the file is open */
    BOOLEAN defined;    /* set true if headP and chanP are set */
    BOOLEAN updateHead; /* set TRUE to force update on close */
    BOOLEAN bReadOnly;  /* if TRUE, no writes allowed */
#if defined(macintosh) || defined(_MAC)
    int     refNum ;    /* file refnum - handle means something else on Mac */
#else
  #ifdef LLIO
    int     handle;     /* file handle */
  #else
    FILE    *handle;    /* stream identifier */
  #endif
#endif
#ifdef USEHANDLES       /* Items required for memory handles */
    THandle headH ;     /* handle to area for file head */
    THandle chanH ;     /* handle to area for channels */
    THandle bufferH ;   /* handle to area to read/write via */
    THandle speedH;     /* handle to area for speed ptrs */
#endif
    TpFileHead headP;   /* pointer to area for file head */
    TpChannel chanP;    /* pointer to area for channels */
    BOOLEAN buffSet;    /* True if buffer available */
    TpDataBlock bufferP;/* pointer to area to read/write via */
    TpSpeedPtr speedP;  /* pointer to area for speed ptrs */
    WORD lastchanRead;  /* Last chan we did a readblock from */
} TSonFile;


extern TSonFile _near files[MAXFILES];  /* the files */
extern TpDataBlock _near workP; /* points at DISKBLOCK bytes work area */

#if (defined(_IS_MSDOS_) || defined(_IS_WINDOWS_)) && !defined(_MAC)
#pragma pack()
#endif

#ifdef __cplusplus
extern "C" {
#endif


/*
** Now declarations of the functions defined in the code
*/
SONAPI(short) SONRead(short fh, void FAR * buffer, WORD bytes, long offset);
SONAPI(short) SONWrite(short fh, void FAR * buffer, WORD bytes, long offset);
SONAPI(short) SONGetBlock(short fh, long offset);
SONAPI(long) SONGetPred(short fh, long offset) ;
SONAPI(long) SONGetSucc(short fh, long offset);
SONAPI(short) SONSetSucc(short fh, long offset, long succOffs);
SONAPI(long) SONFindBlock(short fh, WORD channel,TSTime sTime,TSTime eTime);
SONAPI(short) SONReadBlock(short fh, WORD channel, long position);
SONAPI(short) SONWriteBlock(short fh );
SONAPI(TpChannel) SONChanPnt(short fh, WORD chan);
SONAPI(TSTime) SONIntlChanMaxTime(short fh, WORD chan);
SONAPI(TSTime) SONIntlMaxTime(short fh);
SONAPI(long) SONUpdateMaxTimes(short fh);
SONAPI(void) SONExtendMaxTime(short fh, long time);
SONAPI(long) SONGetFirstData(short fh);

#if defined(macintosh) || defined(_MAC)
SONAPI(int) SONFileHandle(short fh);
#else
    #ifdef LLIO
    SONAPI(int) SONFileHandle(short fh);
    #else
    SONAPI(FILE *) SONFileHandle(short fh);
    #endif
#endif

#ifdef __cplusplus
}
#endif

#endif



