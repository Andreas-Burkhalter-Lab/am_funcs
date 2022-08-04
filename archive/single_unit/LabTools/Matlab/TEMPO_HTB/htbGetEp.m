function epoch = htbGetEp(htb, nEpoch);
%htbGetEp - Read an epoch from file
%
% IN
%   htb             htb structure returned by htbGetHd()
%   nEpoch          Epoch to read (=1, ..., htb.sweep)
%                   For SUM, XSUM, USUM and ESUM data, nEpoch must be 1
%                   For APP, XAPP, UAPP and EAPP data, nEpoch is [1..htb.sweep]
%
% OUT
%   epoch           An N x M matrix of data
%                   The columns are channels (htb.nchannels columns)
%                   The rows are time series of epoch data (htb.period rows)
%                   If htb.sweep==0, a null matrix is returned

epoch = [];                       % No epoch data yet

if (htb.sweep <= 0)
    return;                       % No epoch data
    end

switch htb.func
    case 0                        %  HTYPE_SUM;   SUM for 8 bit ANALOG data
        format = 'int16';
        bytesPerSample = 2;
        maxEpoch = 1;
    case 1                        %  HTYPE_APP;   APP for 8 bit ANALOG data
        format = 'int8';
        bytesPerSample = 1;
        maxEpoch = htb.sweep;
    case 2                        %  HTYPE_USUM;  SUM for COUNTER data
        format = 'uint16';
        bytesPerSample = 2;
        maxEpoch = 1;
    case 3                        %  HTYPE_UAPP;  APP for COUNTER data
        format = 'uint16';
        bytesPerSample = 2;
        maxEpoch = htb.sweep;
    case 4                        %  HTYPE_ESUM;  SUM for EVENT data
        format = 'uint16';
        bytesPerSample = 2;
        maxEpoch = 1;
    case 5                        %  HTYPE_EAPP;  APP for EVENT data
        format = 'uint16';
        bytesPerSample = 2;
        maxEpoch = htb.sweep;
    case 6                        %  HTYPE_XSUM;  XSUM for 12 bit analog data
        format = 'int16';
        bytesPerSample = 2;
        maxEpoch = 1;
    case 7                        %  HTYPE_XAPP;  XAPP for 12 bit analog data
        format = 'int16';
        bytesPerSample = 2;
        maxEpoch = htb.sweep;
    otherwise 
        return;                   % Unknown database format!
    end

if ((nEpoch < 1) | (nEpoch > maxEpoch))
    return;                       % Invalid epoch number
    end

nCols = htb.nchannels;
nRows = htb.period;

% COMPUTE FILE OFFSET TO BEGINNING OF REQUESTED EPOCH

fileoffset = htb.fileoffset;            % htbHeader
fileoffset = fileoffset + 512;          % 512=sizeof(htbHeader)
fileoffset = fileoffset + (bytesPerSample * nCols * nRows) * (nEpoch - 1);

status = fseek(htb.fid, fileoffset, 'bof');   % position to beginning of epoch
if (status == -1)
    return;                             % Error seeking file
    end

% READ EPOCH DATA INTO A MATRIX

[epoch, count] = fread(htb.fid, [nCols nRows], format);

% THE fread() FUNCTION READ THE DATA SO THAT THE ROWS ARE
% THE TIME AND THE COLUMNS ARE THE CHANNELS.  THE HTB.EXE
% UTILITY PROGRAM DISPLAYS THE DATA WHERE THE COLUMNS ARE
% CHANNELS AND THE ROWS ARE TIME.  SO WE TRANSPOSE THE MATRIX.
%
% TRANSPOSE ARRAY SO THAT THE COLUMNS ARE THE CHANNELS AND
% THE ROWS ARE TIME.

epoch = epoch';

return;

 