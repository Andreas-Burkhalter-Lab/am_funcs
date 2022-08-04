unit SON32Imp;		{ SON32.DLL Import-Unit }

interface

uses Windows;

const
  SONMAXCHANS = 32;

  NO_FILE = -1;
  NO_DOS_FILE = -2;
  NO_PATH = -3;
  NO_HANDLES = -4;
  NO_ACCESS = -5;
  BAD_HANDLE = -6;
  MEMORY_ZAP = -7;
  OUT_OF_MEMORY = -8;
  INVALID_DRIVE = -15;
  OUT_OF_HANDLES = -16;
  FILE_ALREADY_OPEN = -600;
  BAD_READ = -17;
  BAD_WRITE = -18;

  NO_CHANNEL = -9;
  CHANNEL_USED = -10;
  CHANNEL_UNUSED = -11;
  PAST_EOF = -12;
  WRONG_FILE = -13;
  NO_EXTRA = -14;
  CORRUPT_FILE = -19;
  PAST_SOF = -20;

  FastWrite:short = 0;
  NormalWrite:short = 1;

type
  TDataKind = (ChanOff,      { = 0}
               Adc,          { = 1}
               EventFall,    { = 2}
               EventRise,    { = 3}
               EventBoth,    { = 4}
               Marker,       { = 5}
               AdcMark,      { = 6, this is marker plus Adc data attached }
               RealMark,     { = 7, Marker with real numbers attached }
               TextMark);    { = 8, Marker with text attached }

  TSTime=longint;
  TAdc=short;
  TMarkBytes=array[0..3] of byte;
  TMarker=record
    mark:TSTime;
    mvals:TMarkBytes;
  end;
  TAdcMark=record
    m:TMarker;
    a:array[0..7999] of TAdc;
  end;
  TRealMark=record
    m:TMarker;
    r:array[0..79] of single;
  end;
  TTextMark=record
    m:TMarker;
    t:array[0..249] of char;
  end;
  TStr=array[0..255] of char;

  TpAdc=^TAdc;
  TpSTime=^TSTime;
  TpMarker=^TMarker;
  TpAdcMark=^TAdcMark;
  TpStr=^TStr;
  TpWord=^Word;
  TpBool=^Boolean;
  TpFloat=^Single;
  TpShort=^Short;

  TFilterElt=char;
  TLayerMask=array[0..31] of TFilterElt;

  TFilterMask=record
    cAllSet:array[0..3] of char;
    aMask  :array[0..3] of TLayerMask;
  end;

  TpFilterMask=^TFilterMask;

procedure SONInitFiles; stdcall;
function SONOpenOldFile(FileName:string):short; stdcall;
function SONOpenNewFile(FileName:string;fMode:short;extra:word):short; stdcall;
function SONCloseFile(fh:short):short; stdcall;
function SONSetBuffSpace(fh:short):short; stdcall;
procedure SONSetFileClock(fh:short;usPerTime:word;TimePerADC:word); stdcall;
function SONSetADCChan(fh:short; Chan:word; adcChan,dvd,buffSz:short; Comment,title:TpStr;
                       ideal,Scale,Offset:single; Units:TpStr):short; stdcall;
function SONSetEventChan(fh:short; Chan:word; evtChan:short; buffSz:short; Comment,title:TpStr;
                         ideal:single; evtKind:TDataKind):short; stdcall;
procedure SONSetFileComment(fh:short;which:word;Comment:string); stdcall;
procedure SONGetFileComment(fh:short;which:word;Comment:TpStr;Max:short); stdcall;
procedure SONSetChanComment(fh:short;Chan:word;Comment:string); stdcall;
procedure SONGetChanComment(fh:short;Chan:word;Comment:TpStr;Max:short); stdcall;
procedure SONSetChanTitle(fh:short;Chan:word;Title:string); stdcall;
procedure SONGetChanTitle(fh:short;Chan:word;Title:TpStr); stdcall;
procedure SONGetIdealLimits(fh:short;Chan:word;ideal,min,max:TpFloat); stdcall;
function SONGetusPerTime(fh:short):word; stdcall;
function SONGetTimePerADC(fh:short):word; stdcall;
procedure SONGetADCInfo(fh:short;Chan:word;Scale:TpFloat;Offset:TpFloat;Units:TpStr;
                        Points:TpWORD;preTrig:TpShort); stdcall;
function SONWriteEventBlock(fh:short;Chan:word;buffer:TpSTime;
                            count:longint):short; stdcall;
function SONWriteMarkBlock(fh:short;Chan:word;buffer:TpMarker;
                           count:longint):short; stdcall;
function SONWriteADCBlock(fh:short;Chan:word;buffer:TpAdc;
                          count:longint;sTime:TSTime):TSTime; stdcall;
function SONGetEventData(fh:short;Chan:word;evtDataP:TpSTime;max:word;sTime,eTime:TSTime;
                         levLowP:TpBOOL;pFiltMask:TpFilterMask):short; stdcall;
function SONGetMarkData(fh:short;Chan:word;markP:TpMarker;max:word;sTime,eTime:TSTime;
                        pFiltMask:TpFilterMask):short; stdcall;
function SONGetADCData(fh:short;Chan:word;adcDataP:TpAdc;max:word;sTime,eTime:TSTime;

                       pbTime:TpSTime;pFiltMask:TpFilterMask):short; stdcall;
function SONChanKind(fh:short;Chan:word):TDataKind; stdcall;
function SONChanDivide(fh:short;Chan:word):TSTime; stdcall;
function SONMaxTime(fh:short):TSTime; stdcall;


implementation

procedure SONInitFiles; stdcall; external 'SON32.DLL';
function SONOpenOldFile(FileName:string):short; stdcall; external 'SON32.DLL';
function SONOpenNewFile(FileName:string;fMode:short;extra:word):short; stdcall; external 'SON32.DLL';
function SONCloseFile(fh:short):short; stdcall; external 'SON32.DLL';
function SONSetBuffSpace(fh:short):short; stdcall; external 'SON32.DLL';
procedure SONSetFileClock(fh:short;usPerTime:word;TimePerADC:word); stdcall; external 'SON32.DLL';
function SONSetADCChan(fh:short; Chan:word; adcChan,dvd,buffSz:short; Comment,title:TpStr;
                       ideal,Scale,Offset:single; Units:TpStr):short; stdcall; external 'SON32.DLL';
function SONSetEventChan(fh:short; Chan:word; evtChan:short; buffSz:short; Comment,title:TpStr;
                         ideal:single; evtKind:TDataKind):short; stdcall; external 'SON32.DLL';
procedure SONSetFileComment(fh:short;which:word;Comment:string); stdcall; external 'SON32.DLL';
procedure SONGetFileComment(fh:short;which:word;Comment:TpStr;Max:short); stdcall; external 'SON32.DLL';
procedure SONSetChanComment(fh:short;Chan:word;Comment:string); stdcall; external 'SON32.DLL';
procedure SONGetChanComment(fh:short;Chan:word;Comment:TpStr;Max:short); stdcall; external 'SON32.DLL';
procedure SONSetChanTitle(fh:short;Chan:word;Title:string); stdcall; external 'SON32.DLL';
procedure SONGetChanTitle(fh:short;Chan:word;Title:TpStr); stdcall; external 'SON32.DLL';
procedure SONGetIdealLimits(fh:short;Chan:word;ideal,min,max:TpFloat); stdcall; external 'SON32.DLL';
function SONGetusPerTime(fh:short):word; stdcall; external 'SON32.DLL';
function SONGetTimePerADC(fh:short):word; stdcall; external 'SON32.DLL';
procedure SONGetADCInfo(fh:short;Chan:word;Scale:TpFloat;Offset:TpFloat;Units:TpStr;
                        Points:TpWORD;preTrig:TpShort); stdcall; external 'SON32.DLL';
function SONWriteEventBlock(fh:short;Chan:word;buffer:TpSTime;
                            count:longint):short; stdcall; external 'SON32.DLL';
function SONWriteMarkBlock(fh:short;Chan:word;buffer:TpMarker;
                           count:longint):short; stdcall; external 'SON32.DLL';
function SONWriteADCBlock(fh:short;Chan:word;buffer:TpAdc;
                          count:longint;sTime:TSTime):TSTime; stdcall; external 'SON32.DLL';
function SONGetEventData(fh:short;Chan:word;evtDataP:TpSTime;max:word;sTime,eTime:TSTime;
                         levLowP:TpBOOL;pFiltMask:TpFilterMask):short; stdcall; external 'SON32.DLL';
function SONGetMarkData(fh:short;Chan:word;markP:TpMarker;max:word;sTime,eTime:TSTime;
                        pFiltMask:TpFilterMask):short; stdcall; external 'SON32.DLL';
function SONGetADCData(fh:short;Chan:word;adcDataP:TpAdc;max:word;sTime,eTime:TSTime;
                       pbTime:TpSTime;pFiltMask:TpFilterMask):short; stdcall; external 'SON32.DLL';
function SONChanKind(fh:short;Chan:word):TDataKind; stdcall; external 'SON32.DLL';
function SONChanDivide(fh:short;Chan:word):TSTime; stdcall; external 'SON32.DLL';
function SONMaxTime(fh:short):TSTime; stdcall; external 'SON32.DLL';

end.

