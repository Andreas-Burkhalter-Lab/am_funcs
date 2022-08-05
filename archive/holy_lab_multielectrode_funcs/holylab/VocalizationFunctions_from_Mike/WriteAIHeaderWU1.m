function filePosition = WriteAIHeaderWU1(fid,h,headerType)
% WRITEAIHEADERWU1: Write the header for legacy raw waveform files
%
% This function is called to write the header to many different file types.
% Each file type's header has univeral information (AI header type) as well as 
% information that is specific to that particular file type.  The universal 
% information is version 3. The specific information is version 1.
%
% Syntax:
%    filePosition = WriteAIHeaderWU1(fid,h,headerType)
% where
%    fid is the numeric file identifier
%    h is a structure that contains the information that is to be written to the header
%    headerType is a structure that specifies which type of header is to be written
% and
%    filePosition is a structure that returns the position in the header for various
%      elements that need to be updated
%
% See also: READAIHEADERWU1

magicNum = 992732356;
fwrite(fid,magicNum,'int32');
headersize = 0;                        % Update this later
fwrite(fid,headersize,'uint32');
type = 0;
% Set individual bits to true if a particular header type is to be written
if (isfield(headerType,'AI') & headerType.AI == 1)
    type = bitset(type,1);
end

if (isfield(headerType,'Sng') & headerType.Sng == 1)
    type = bitset(type,4);
end

if (isfield(headerType,'Detect') & headerType.Detect == 1)
    type = bitset(type,5);
end    
% Put checks for other header types here
fwrite(fid,type,'int32');
if (bitget(type,1))
        % Write universal AI header type
        version = 3;
        fwrite(fid,version,'int16');        % Header version
        filePosition.nscans = ftell(fid);
        fwrite(fid,h.nscans,'uint32');
        numch = prod(size(h.channels));
        fwrite(fid,numch, 'int32');
        fwrite(fid,h.channels,'int16');
        fwrite(fid,h.scanrate,'float32');
        fwrite(fid,h.scalemult,'float32');
        fwrite(fid,h.scaleoff,'float32');
        fwrite(fid,h.voltageMin,'int16');
        fwrite(fid,h.voltageMax,'int16');
        WriteLVString(fid,h.date);
        WriteLVString(fid,h.time);
        WriteLVString(fid,h.usrhdr);
end
if (bitget(type,4))
    % Write sonogram header type
    version = 1;
    fwrite(fid,version,'int16');
    fwrite(fid,h.nfreq + 1,'int32');
    filePosition.columnTotal = ftell(fid);
    fwrite(fid,h.columnTotal,'int32');
    fwrite(fid,h.threshold,'float32');
    filePosition.nblocks = ftell(fid);
    fwrite(fid,h.nblocks,'int32');
    filePosition.tacq = ftell(fid);
    fwrite(fid,h.tacq,'float32');
    fwrite(fid,h.freqMin,'int32');
    fwrite(fid,h.freqMax,'int32');
end   
if (bitget(type,5))
    % Write proximity detection header type
    version = 1;
    fwrite(fid,version,'int16');
    filePosition.numTransitions = ftell(fid);
    fwrite(fid,h.numTransitions,'int16');
end

% Update header size
fcur = ftell(fid);
fseek(fid,4,-1);
fwrite(fid,fcur,'uint32');
fseek(fid,fcur,-1);
