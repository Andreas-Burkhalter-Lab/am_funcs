function err = htbWrite(fid, htb, epochData)
%htbWrite - Append a database to end of an HTB file
%
% The HTB file must be opened for write with htbOpenW().
% If the HTB file doesn't exist, it is created.
% The header and database are appended to the end of the HTB file.
% Other databases in the file are not changed.
% No validity checks are made on the header or data.
%
% SYNOPSIS
%   err = htbWrite(fid, htb, epochData);
%
% IN
%   fid             fid returned by htbOpenw()
%   htb             htb structure
%   epochData       an N x M matrix where M=# of channels
%                       Append databases N = htb.nsweeps * htb.period
%                       Sum databases    N = htb.nsweeps
%
% OUT
%   err             0 if successful, non-zero if error

% Note 1: We are careful to write bytes sequentially because the file
% was opened with 'a+' which is "append" attribute by htbOpenW().
% Note 2: We do not worry here about he endian problem.  We assume
% that this code is running on a computer that is native little-endian
% as the IBM PCs are.

err = AppendHeader(fid, htb);
if (err == 0)
    err = AppendData(fid, htb, epochData);
    end
return;

%---------------------------------------------------------------------
function err = AppendHeader(fid, htb)

PutSTRING(fid, htb.date, 26);
PutLONG(fid, htb.ldate);
PutSTRING(fid, htb.cfg_file, 14);
PutSTRING(fid, htb.pro_file, 14);
PutZeroBytes(fid, 52);
PutULONG(fid, htb.speed);
PutULONG(fid, htb.alloc);

PutLONG(fid, htb.offset);
PutULONG(fid, htb.period);
PutULONG(fid, htb.extension);
PutUSHORT(fid, htb.skip);
PutUSHORT(fid, htb.first_channel);
PutUSHORT(fid, htb.nchannels);
PutUSHORT(fid, htb.sweep_limit);
PutULONG(fid, htb.cancel_override);
PutUCHAR(fid, htb.func);

PutZeroBytes(fid, 1);
PutUSHORT(fid, htb.tag);

PutUSHORT(fid, htb.npages);
PutULONG(fid, htb.nsamples);
PutUSHORT(fid, htb.samples_per_page);
PutUSHORT(fid, htb.sweep);
PutUSHORT(fid, htb.next_page);
PutUSHORT(fid, htb.next_off);

PutSTRING(fid, htb.title, 80);
PutULONG(fid, htb.speed_units);
PutZeroBytes(fid, 268);

err = 0;
return;

%---------------------------------------------------------------------
function err = AppendData(fid, htb, epochData)

epochData = epochData';           % fwrite writes columns first

switch htb.func
    case 0                        %  HTYPE_SUM;   SUM for 8 bit ANALOG data
        format = 'int16';
    case 1                        %  HTYPE_APP;   APP for 8 bit ANALOG data
        format = 'int8';
    case 2                        %  HTYPE_USUM;  SUM for COUNTER data
        format = 'uint16';
    case 3                        %  HTYPE_UAPP;  APP for COUNTER data
        format = 'uint16';
    case 4                        %  HTYPE_ESUM;  SUM for EVENT data
        format = 'uint16';
    case 5                        %  HTYPE_EAPP;  APP for EVENT data
        format = 'uint16';
    case 6                        %  HTYPE_XSUM;  XSUM for 12 bit analog data
        format = 'int16';
    case 7                        %  HTYPE_XAPP;  XAPP for 12 bit analog data
        format = 'int16';
    otherwise 
        printf('AppendData - Unknown htb.func = %d', htb.func);
        return;                   % Unknown database format!
    end

count = fwrite(fid, epochData, format);
err = 0;
return;


%-----------------------------------------------------------------------------
function err = PutZeroBytes(fid, nBytes)
%PutZeroBytes - write a number of 0 bytes

count = fwrite(fid, linspace(0,0,nBytes), 'uint8');
if (count ~= nBytes)
    err = -1;
    end;

return;

%-----------------------------------------------------------------------------
function err = PutUCHAR(fid, cData)
%PutUCHAR - Write unsigned char to HTB file

err = 0;
count = fwrite(fid, cData, 'uint8');

if (count ~= 1)
    err = -1;
    end;
return;

%-----------------------------------------------------------------------------
function err = PutUSHORT(fid, nUSHORT)
%PutUSHORT - Write unsigned short to HTB file

err = 0;
count = fwrite(fid, nUSHORT, 'uint16');

if (count ~= 1)
    err = -1;
    end;
return;

%-----------------------------------------------------------------------------
function err = PutULONG(fid, nULONG)
%PutULONG - Write unsigned long to HTB file

err = 0;
count = fwrite(fid, nULONG, 'uint32');

if (count ~= 1)
    err = -1;
    end;
return;

%-----------------------------------------------------------------------------
function err = PutLONG(fid, nLONG)
%PutLONG - Write long to HTB file

err = 0;
count = fwrite(fid, nLONG, 'int32');

if (count ~= 1)
    err = -1;
    end;
return;

%-----------------------------------------------------------------------------
function err = PutSTRING(fid, pSTRING, nbytes)
%PutSTRING - write string to HTB file
% 
% We convert the string to an array of bytes.
% Then we right pad with 0 bytes, if necessary.
% Then we trim the array to exactly nbytes and write it out.

err = 0;
t = double(pSTRING);                    % Convert string to array of numbers
t(nbytes+1) = 0;                        % Right pad with 0s so that it is at least nbytes+1
t = t(1:nbytes);                        % Trim it to exactly nbytes
count = fwrite(fid, t, 'int8');         % Write it out as an array of bytes

if (count ~= nbytes)
    err = -1;
    end;
return;
