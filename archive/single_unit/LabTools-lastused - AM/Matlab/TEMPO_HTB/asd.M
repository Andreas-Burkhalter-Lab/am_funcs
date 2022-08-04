function htb = htbGetHd(fid, db)
%htbGetHd - Read HTB database db's header
%
% IN
%   fid             fid returned by htbOpen() or htbOpenw()
%   db              A database number in file =1,...,htbCount(fid)
%
% OUT
%   htb             An HTB header structure
%                   htb.alloc == 0 if an error occured
%
% The returned value is an HTB structure with the following fields:
%
%    char    date[26];                   time/date saved (formatted by ctime)
%    long    ldate;                      time/date (from time) equivalent to date
%    char    cfg_file[14];               Null terminated configuration file name
%    char    pro_file[14];               Null terminated protocol file name
%    ULONG   speed;                      Speed of acquisition
%    ULONG   alloc;                      # of bytes occupied by database in file (including header)
%
%    long    offset;                     Offset from trigger point 
%    ULONG   period;                     length of histogram 
%    ULONG   extension;                  cancel vunerable time after period 
%    USHORT  skip;                       # of ring buffer bins between histogram bins (0..)
%    USHORT  first_channel;              First channel to accumulate (0, 1, .. ) 
%    USHORT  nchannels;                  # of channels in a channel set (1, 2, .. ) 
%    USHORT  sweep_limit;                Maximum number of functions 
%    ULONG   cancel_override;            Cancel override bits 
%    UCHAR   func;                       Function to perform (0,1,...)
%    USHORT  tag;                        Class to which this histogram is a member 
%
%    // THE FOLLOWING ARE COMPUTED FROM THE ABOVE VALUES
%    // These values must NOT be changed except by the KERNEL
%        
%    USHORT  npages;                     # of 16kb EMM pages allocated 
%    ULONG   nsamples;                   Total number of channel sets 
%    USHORT  samples_per_page;           # of channel sets per page 
%    USHORT  sweep;                      Current number of sweeps 
%    USHORT  next_page, next_off;        Next location for new data 
%
%    char    title[80];                  H Table's title
%    ULONG   speed_units;                Units of speed

htb = htbInit;                          % Get an empty HTB structure

offset = GetOffset(fid, db);            % Offset to HTB header in file
if (offset == -1)
    return;
end;
    
[htb.alloc, count] = GetULONG(fid, offset + 114);  % 114 = htb.alloc

if ((count ~= 1) | (htb.alloc <= 0))
    return;
    end;

[htb.date, count]           = GetString(fid, offset + 0, 26);
[htb.ldate, count]          = GetLONG(fid, offset + 26);
[htb.cfg_file, count]       = GetString(fid, offset + 30, 14);
[htb.pro_file, count]       = GetString(fid, offset + 44, 14);
[htb.speed, count]          = GetULONG(fid, offset + 110);
[htb.offset, count]         = GetLONG(fid, offset + 118);
[htb.period, count]         = GetULONG(fid, offset + 122);
[htb.extension, count]      = GetULONG(fid, offset + 126);
[htb.skip, count]           = GetUSHORT(fid, offset + 130);
[htb.first_channel, count]  = GetUSHORT(fid, offset + 132);
[htb.nchannels, count]      = GetUSHORT(fid, offset + 134);
[htb.sweep_limit, count]    = GetUSHORT(fid, offset + 136);
[htb.cancel_override, count]= GetULONG(fid, offset + 138);
[htb.func, count]           = GetUCHAR(fid, offset + 142);
[htb.tag, count]            = GetUSHORT(fid, offset + 144);
[htb.npages, count]         = GetUSHORT(fid, offset + 146);
[htb.nsamples, count]       = GetULONG(fid, offset + 148);
[htb.samples_per_page, count]= GetUSHORT(fid, offset + 152);
[htb.sweep, count]          = GetUSHORT(fid, offset + 154);
[htb.next_page, count]      = GetUSHORT(fid, offset + 156);
[htb.next_off, count]       = GetUSHORT(fid, offset + 158);
[htb.title, count]          = GetString(fid, offset + 160, 80);
[htb.speed_units, count]    = GetULONG(fid, offset + 240);

htb.fileoffset = offset;    % Remember file offset for later use
htb.fid = fid;              % Remember file id for later use

return;


%-----------------------------------------------------------------------------
function offset = GetOffset(fid, db)
%GetOffset - Returns file offset of HTB header for db

offset = 0;                             % Offset of first HTB structure in file
d = 1;
while (d < db)
    [nbytes, count] = GetULONG(fid, offset + 114); % 114=offset to htb.alloc field
    
    if ((count ~= 1) | (nbytes <= 0))
        offset = -1;
        return;
        end;
    
    offset = offset + nbytes;
    d = d + 1;        
    end;
return;


%-----------------------------------------------------------------------------
function [i, count] = GetUCHAR(fid, offset)
%GetUCHAR - Read unsigned char from HTB file

i = 0;
count = -1;

status = fseek(fid, offset, 'bof');
if (status == -1)
    offset = -1;
    return;
    end;

i = 0;    
[i, count] = fread(fid, 1, 'uint8');  % Read htb.alloc

if (count ~= 1)
    i = 0;
    return;
    end;
return;

%-----------------------------------------------------------------------------
function [i, count] = GetUSHORT(fid, offset)
%GetUSHORT - Read unsigned short from HTB file

i = 0;
count = -1;

status = fseek(fid, offset, 'bof');
if (status == -1)
    offset = -1;
    return;
    end

i = 0;    
[i, count] = fread(fid, 1, 'uint16');  % Read htb.alloc

if (count ~= 1)
    i = 0;
    return;
    end;
return;

%-----------------------------------------------------------------------------
function [i, count] = GetULONG(fid, offset)
%GetULONG - Read unsigned long from HTB file

i = 0;
count = -1;

status = fseek(fid, offset, 'bof');
if (status == -1)
    offset = -1;
    return;
    end

i = 0;    
[i, count] = fread(fid, 1, 'uint32');  % Read htb.alloc

if (count ~= 1)
    i = 0;
    return;
    end;
return;

%-----------------------------------------------------------------------------
function [i, count] = GetLONG(fid, offset)
%GetLONG - Read long from HTB file

i = 0;
count = -1;

status = fseek(fid, offset, 'bof');
if (status == -1)
    offset = -1;
    return;
    end

i = 0;    
[i, count] = fread(fid, 1, 'int32');  % Read htb.alloc

if (count ~= 1)
    i = 0;
    return;
    end;
return;

%-----------------------------------------------------------------------------
function [i, count] = GetString(fid, offset, nbytes)
%GetString - Read long from HTB file

i = '';
count = -1;

status = fseek(fid, offset, 'bof');
if (status == -1)
    offset = -1;
    return;
    end

i = '';
[i, count] = fread(fid, [1 nbytes], 'char');  % Read htb.alloc
i = char(i);
i = deblank(i);
if (count ~= nbytes)
    return;
    end;
return;

