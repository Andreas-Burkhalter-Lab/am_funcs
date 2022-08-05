function [header,headersize] = ReadAIHeaderHarvard(file)
% ReadAIHeaderHarvard: Read the header from legacy raw waveform files
% [header,headersize] = ReadAIHeaderHarvard(file)
% AI files are written either by (version 1) the
%        LabVIEW VI 'create AI binary header.vi,' or
%        (version 2) by C code.
% Other header creators are derived from this one, and
% start with it.
% If "file" is a string, it's treated as a filename
% If "file" is numeric, it's treated as a file identifier
if (ischar(file))
        [fid,message] = fopen(file,'r');
        if (fid < 1)
                error(message);
        end
elseif (isnumeric(file))
        fid = file;
else
        error(['Do not recognize input ',file]);
end
headersize = fread(fid,1,'uint32');
header.headersize = headersize;
header.type = fread(fid,1,'int16');
if (~(header.type == 1 | header.type == 4 | header.type == 8))
        error('File is not AI type');
end
header.version = fread(fid,1,'int16');
if (header.version == 1)
        header.chstr = readLVstring(fid);
        chconfiglen = fread(fid,1,'uint32');
        header.numch = fread(fid,1,'int32');
        for i=1:header.numch
                header.chinfo(i) = readChConfig(fid);
        end
        header.scanrate = fread(fid,1,'float32');
        header.chclock = fread(fid,1,'float32');
        header.date = readLVstring(fid);
        header.time = readLVstring(fid);
        header.usrhdr = readLVstring(fid);
        header.acqtime = fread(fid,1,'float32');
        header.nscans = fread(fid,1,'uint32');
        header.scalemult = header.chinfo(1).scalemult;
        header.scaleoff = header.chinfo(1).scaleoff;
        header.voltageMin = -round(4095*header.scalemult/2);
        header.voltageMax = -header.voltageMin;
        header.channels = str2num(header.chstr');
        if (isempty(header.channels))
                for i = 1:header.numch
                        header.channels(i) = str2num(header.chinfo(i).channel);
                end
        end
        if (size(header.channels,2) == 1)
                header.channels = header.channels';
        end
        header.minSample = -2048;
        header.maxSample = 2047;
        header.flag = 'Harvard';
else
        header.nscans = fread(fid,1,'uint32');
        header.numch = fread(fid,1,'int32');
        header.channels = fread(fid,header.numch,'int16');
        if (size(header.channels,2) == 1)
                header.channels = header.channels';
        end
        header.scanrate = fread(fid,1,'float32');
        header.scalemult = fread(fid,1,'float32');
        header.scaleoff = fread(fid,1,'float32');
        header.voltageMin = -round(4095*header.scalemult/2);
        header.voltageMax = -header.voltageMin;
        header.date = ReadLVString(fid);
        header.time = ReadLVString(fid);
        header.usrhdr = ReadLVString(fid);
        header.minSample = -2048;
        header.maxSample = 2047;
        header.flag = 'Harvard';
end
if (ischar(file))
        fclose(fid);
end
