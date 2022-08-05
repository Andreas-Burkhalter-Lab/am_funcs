function scan = msload_mzXML(filename,options)
% msload_mzXML: load mass spectrometry data from mzXML format
% Syntax:
%   scan = msload_mzXML(filename)
%   scan = msload_mzXML(filename,options)
% where
%   filename is a string containing the full name of the file
%   options is a structure which may have the following fields:
%     ms_level (default 'all'): if a numeric value is supplied, only scans
%       with the given level are loaded (MS = 1, MS/MS = 2, MS/MS/MS = 3,
%       etc.)
%     scan_time_range (default [0 Inf]): the time range of scans to load,
%       in seconds
%     scan_number_range: if supplied, loads only scans within these indices
% and
%   scan is a structure array, each element containing data for a single
%     scan.

% Copyright 2009 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  options = default(options,'ms_level','all','scan_time_range',[0 Inf]);
  [fid,msg] = fopen(filename);
  if (fid < 0)
    error(msg);
  end
  
  % Go to the end of the file and extract the location of the index
  tail_offset = 200;
  status = fseek(fid,-tail_offset,'eof');
  if (status < 0)
    error('File was not long enough to be mzXML');
  end
  tailstr = fread(fid,[1 inf],'*char');

  indx = passtoken(tailstr,'<indexOffset>');
  if isempty(indx)
    error('Can''t find the index');
  end
  index_position = sscanf(tailstr(indx:end),'%d');
  
  % Go to the index and read it in
  fseek(fid,index_position,'bof');
  indexstr = fread(fid,[1 inf],'*char');
  scan_offset_indx_start = regexp(indexstr,'<offset id="[0-9]+" >','end')+1;
  scan_offset_indx_end = regexp(indexstr,'</offset>','start')-1;
  n_scans = length(scan_offset_indx_start);
  scan_offset = zeros(1,n_scans);
  for i = 1:n_scans
    scan_offset(i) = str2double(indexstr(scan_offset_indx_start(i):scan_offset_indx_end(i)));
  end
  % Add in the beginning of the index as a sentinel
  scan_offset = [scan_offset index_position];
  
  % Read the requested scans
  scanrange = [1 n_scans];
  if isfield(options,'scan_number_range')
    scanrange = options.scan_number_range;
  end
  scans_to_read = scanrange(1):scanrange(end);
  n_scans_to_read = length(scans_to_read);
  scan_len = diff(scan_offset);
  scanC = cell(1,n_scans_to_read);
  cIndex = 1;
  for i = 1:n_scans_to_read
    thisScan = scans_to_read(i);
    fseek(fid,scan_offset(thisScan),'bof');
    scanstr = fread(fid,[1 scan_len(thisScan)],'*char');
    tmp = parse_scanstr(scanstr,options);
    if ~isempty(tmp)
      scanC{cIndex} = tmp;%   str = uint8(peaksdatac);
%   n_peaks = n_points;
%   save test_in str n_peaks

      cIndex = cIndex+1;
    end
  end
  scan = cat(2,scanC{1:cIndex-1});
  fclose(fid);
  
  % Look for corrupted scans
  n_scans = length(scan);
  goodscan = true(1,n_scans);
  for i = 1:n_scans
    if any(scan(i).mz < 0.99*scan(i).min_mz | scan(i).mz > 1.01*scan(i).max_mz)
      goodscan(i) = false;
    end
  end
  scan = scan(goodscan);
  if any(~goodscan)
    warning('mzXML:corruptedscan','%d corrupted scans detected',sum(~goodscan));
  end
end


function s = kv(str)
  mat = regexp(str,'\w*=\S*','match');
  for i = 1:length(mat)
    indx = regexp(mat{i},'=');
    s.(mat{i}(1:indx-1)) = mat{i}(indx+2:end-1);
  end
end


function indx = passtoken(str,tok)
  indx = strfind(str,tok);
  if isempty(indx)
    return
  end
  indx = indx+length(tok);
end

function scan = parse_scanstr(str,options)
  % Check to see if this scan needs to be parsed
  scan = [];
  if ~isequal(options.ms_level,'all')
    index = regexp(str,'msLevel="','end','once') + 1;
    ms_level = sscanf(str(index:index+10),'%d');
    if (ms_level ~= options.ms_level)
      return;
    end
  end
  index = regexp(str,'retentionTime="','end','once') + 3;
  scan_time = sscanf(str(index:index+10),'%f');
  if (scan_time < options.scan_time_range(1) || scan_time > options.scan_time_range(2))
    return
  end
  
  % Find markers for beginning and end of "blocks" in the header
  header_end = regexp(str,'>','once');
  data_start = regexp(str,'"m/z-int.*" >','once','end')+1;
  [precursor_start,precursor_startmz] = regexp(str(1:data_start),'<precursorMz.*?>','start','end','once');
  
  % Parse key-value pairs
  scanheader = kv(str(1:header_end));
  precursormz = nan;
  precursorI = nan;
  if ~isempty(precursor_start)
    precursor = kv(str(precursor_start:precursor_startmz));
    precursorI = str2double(precursor.precursorIntensity);
    precursormz = sscanf(str((1:20)+precursor_startmz),'%g');
  end
  
  tail_offset = 100;
  scan_tail = str(end-tail_offset:end);
  data_end = regexp(scan_tail,'</peaks>','start');
  if isempty(data_end)
    error('Can''t find end of peaks data');
  end
  data_end = tail_offset-data_end+2;
  peaksdatac = str(data_start:end-data_end);
  
  n_points = str2double(scanheader.peaksCount);
  peaksdata = mzxml_decode(uint8(peaksdatac),n_points);
  
  scan = struct('scan_time',scan_time,...
    'ms_level',str2double(scanheader.msLevel),...
    'min_mz',str2double(scanheader.lowMz),...
    'max_mz',str2double(scanheader.highMz),...
    'totIntensity',str2double(scanheader.totIonCurrent),...
    'precursor_mz',precursormz,...
    'precursor_I',precursorI,...
    'n_points',n_points,...
    'mz',peaksdata(1,:),...
    'intensity',peaksdata(2,:));
end