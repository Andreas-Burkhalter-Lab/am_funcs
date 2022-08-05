function scan = msload_mzML(filename,options)
  % msload_mzML: load mass spectrometry data from mzML format
  % Syntax:
  %   scan = msload_mzML(filename)
  %   scan = msload_mzML(filename,options)
  % where
  %   filename is a string containing the full name of the file
  %   options is a structure which may have the following fields:
  %     ms_level (default 'all'): if a numeric value is supplied, only scans
  %       with the given level are loaded (MS = 1, MS/MS = 2, MS/MS/MS = 3,
  %       etc.)
  %     scan_time_range (default [0 Inf]): the time range of scans to load,
  %       in seconds (must also be included in the scan_number_range below)
  %     scan_number_range: if supplied, loads only scans within these indices
  %       (must also be within the scan_time_range above)
  % and
  %   scan is a structure array, each element containing data for a single
  %     scan.
  
  % Copyright 2012 by Timothy E. Holy
  
  if nargin < 2
    options = struct;
  end
  options = default(options,'ms_level','all','scan_time_range',[0 Inf],'scan_number_range',[1 Inf]);
  
  % We can't use xmlread, because some of these files are huge, and xmlread
  % is very inefficient. It's better to parse the file incrementally.
  ixml = incrementalXML(filename);
  [ixml,cnode] = ixml.next;
  if ~isempty(cnode) || ~strcmp(ixml.open_node{1}.name,'mzML')
    error('Not a mzML file')
  end
  
  %% Parse the 'referenceableParamGroupList'
  str = 'referenceableParamGroupList';
  ixml = ixml.advance(str);
  % Read all referenceableParamGroups
  n = str2double(ixml.current.value('count'));  % # of groups
  msparams = repmat(struct('name','','params',[]),1,n);
  [ixml,cnode] = ixml.next;
  for i = 1:n
    % Read one referenceableParamGroup and store the parameters
    rpg = {};
    str2 = 'referenceableParamGroup';
    while isempty(cnode) || ~strcmp(cnode.name,str2)
      if ~isempty(cnode)
        rpg{end+1} = cnode; %#ok<AGROW>
      end
      [ixml,cnode] = ixml.next;
    end
    msparams(i).name = cnode.value('id');
    msparams(i).params = rpg;
    [ixml,cnode] = ixml.next;
  end
  
  %% Parse the spectra
  scan = [];
  str = 'spectrumList';
  ixml = ixml.advance(str);
  n_scans = str2double(ixml.current.value('count'));  % # of spectra
  options.scan_number_range(2) = min(options.scan_number_range(2),n_scans);
  str2 = 'spectrum';
  for i = 1:options.scan_number_range(1)
    ixml = ixml.advance(str2);
    [ixml,cnode] = ixml.next;
  end
  tic
  tlast = toc;
  first = true;
  n_scans = diff(options.scan_number_range)+1;
  for i = options.scan_number_range(1):options.scan_number_range(2)
    if (toc > tlast + 5)
      tlast = toc;
      if first
        fprintf('%% done: ');
        first = false;
      end
      fprintf('%d...',round(100*i/n_scans))
    end
    mslmzml_read_spectrum;
    if tmp.scan_time < options.scan_time_range(1) || tmp.scan_time > options.scan_time_range(2)
      continue
    end
    if ~isequal(options.ms_level,'all')
      if options.ms_level ~= tmp.ms_level
        continue
      end
    end
    if isempty(scan)
      scan = tmp;
    else
      scan(end+1) = tmp; %#ok<AGROW>
    end
    [ixml,cnode] = ixml.next;
  end
  if ~first
    fprintf('done.\n');
  end
  
  function mslmzml_read_spectrum
    ms_level = [];
    min_mz = [];
    max_mz = [];
    totIntensity = [];
    t = [];
    datastr = '';
    dataname = '';
    compression = '';
    datatype = '';
    datatmp = struct('mz',[],'intensity',[]);
    while isempty(cnode) || ~strcmp(cnode.name,'spectrum')
      if ~isempty(cnode)
        switch cnode.name
          case 'referenceableParamGroupRef'
            % Set parameters
            name = cnode.value('ref');
            pflag = strcmp(name,{msparams.name});
            if sum(pflag) ~= 1
              error('Must be unique match')
            end
            cparams = msparams(pflag).params;
            for pIndex = 1:length(cparams)
              name = cparams{pIndex}.value('name');
              switch name
                case 'ms level'
                  ms_level = str2double(cparams{pIndex}.value('value'));
                case 'scan window lower limit'
                  min_mz = str2double(cparams{pIndex}.value('value'));
                case 'scan window upper limit'
                  max_mz = str2double(cparams{pIndex}.value('value'));
              end
            end
          case 'cvParam'
            name = cnode.value('name');
            switch name
              case 'total ion current'
                totIntensity = str2double(cnode.value('value'));
              case 'elution time'
                t = str2double(cnode.value('value'));
                u = cnode.value('unitName');
                if strcmp(u,'minute')
                  t = t*60;
                end
              case 'm/z array'
                dataname = 'mz';
              case 'intensity array'
                dataname = 'intensity';
              case 'zlib compression'
                compression = 'zlib';
              case '64-bit float'
                datatype = 'double';
            end
          case 'binary'
            datastr = cnode.content;
          case 'binaryDataArray'
            datatmp.(dataname) = xmlbinconvert(datastr,compression,datatype)';
        end
      end
      [ixml,cnode] = ixml.next;
    end
    tmp = struct('id',cnode.value('id'),...
      'scan_time',t,...
      'ms_level',ms_level,...
      'min_mz',min_mz,...
      'max_mz',max_mz,...
      'totIntensity',totIntensity,...
      'precursor_mz',[],...
      'precursor_I',[],...
      'n_points',length(datatmp.mz),...
      'mz',datatmp.mz,...
      'intensity',datatmp.intensity);
  end
end
% 
%     
%   xdoc = xmlread(filename);
% 
%   %% Look for scan parameter lists
%   allparams = xdoc.getElementsByTagName('referenceableParamGroup');
%   n_samples = allparams.getLength/3;  % SpectrumParams, ScanParams, and ScanWindowParams
%   for i = 1:n_samples
%     % Parse SpectrumParams
%     k = 3*(i-1);
%     thisparam = allparams.item(k);
%     msparams(i).spectrumparamsname = char(thisparam.getAttribute('id'));
%     % check that this is profile data and get the ms level
%     allcvp = thisparam.getElementsByTagName('cvParam');
%     thiscvp = allcvp.item(1);
%     msparams(i).ms_level = str2double(char(thiscvp.getAttribute('value')));
%     thiscvp = allcvp.item(2);
%     thisname2 = char(thiscvp.getAttribute('name'));
%     if ~strcmp(thisname2(1:7),'profile')
%       error('Only supports profile data');
%     end
%     % Parse ScanWindowParams
%     k = 3*(i-1)+2;
%     thisparam = allparams.item(k);
%     allcvp = thisparam.getElementsByTagName('cvParam');
%     thiscvp = allcvp.item(0);
%     msparams(i).min_mz = str2double(char(thiscvp.getAttribute('value')));
%     thiscvp = allcvp.item(1);
%     msparams(i).max_mz = str2double(char(thiscvp.getAttribute('value')));
%   end
%   
%   %% Look for "spectra" elements
%   allspectra = xdoc.getElementsByTagName('spectrum');
%   options.scan_number_range(2) = min(options.scan_number_range(2),allspectra.getLength);
%   n_scans = diff(options.scan_number_range)+1;
%   scan = [];
%   for k = 0:n_scans-1
%     thisscannumber = k + options.scan_number_range(1)-1;
%     thisspectrum = allspectra.item(thisscannumber);
%     thisid = char(thisspectrum.getAttribute('id'));
%     % Determine the scan type, so we can look up the parameters
%     allspecparams = thisspectrum.getElementsByTagName('referenceableParamGroupRef');
%     thisspecparams = allspecparams.item(0);
%     ref = char(thisspecparams.getAttribute('ref'));
%     indx = findparams(ref,msparams);  % this is the index into the ms parameters
%     % Check if the ms_level is consistent with the desired type
%     if ~strcmp(options.ms_level,'all') && msparams(indx).ms_level ~= options.ms_level
%       continue; % wrong ms_level, we don't have to read this one
%     end
%     % Get total ion current
%     allcvp = thisspectrum.getElementsByTagName('cvParam');
%     thiscvp = allcvp.item(0);
%     totcurrent = str2double(char(thiscvp.getAttribute('value')));
%     % Get the scan time
%     allscanlist = thisspectrum.getElementsByTagName('scanList');
%     thisscanlist = allscanlist.item(0);
%     allscan = thisscanlist.getElementsByTagName('scan');
%     thisscan = allscan.item(0);
%     allcvp = thisscan.getElementsByTagName('cvParam');
%     thisparam = allcvp.item(0);
%     v = str2double(char(thisparam.getAttribute('value')));
%     u = char(thisparam.getAttribute('unitName'));
%     switch u
%       case 'second'
%         t = v;
%       case 'minute'
%         t = v*60;
%     end
%     % Check that the scan time is within bounds
%     if t < options.scan_time_range(1) || t > options.scan_time_range(2)
%       continue
%     end
%     % Read the data format
%     allbinarylist = thisspectrum.getElementsByTagName('binaryDataArray');
%     if allbinarylist.getLength ~= 2
%       error('Expecting binary list of length 2');
%     end
%     tmp = struct('id',thisid,...
%       'scan_time',t,...
%       'ms_level',msparams(indx).ms_level,...
%       'min_mz',msparams(indx).min_mz,...
%       'max_mz',msparams(indx).max_mz,...
%       'totIntensity',totcurrent,...
%       'precursor_mz',[],...
%       'precursor_I',[],...
%       'n_points',[],...
%       'mz',[],...
%       'intensity',[]);
%     for binaryIndex = 0:1
%       thisbinary = allbinarylist.item(binaryIndex);
%       allcvp = thisbinary.getElementsByTagName('cvParam');
%       fn = '';
%       datatype = '';
%       compression = '';
%       for vecIndex = 0:allcvp.getLength-1
%         thiscvp = allcvp.item(vecIndex);
%         switch char(thiscvp.getAttribute('name'))
%           case 'm/z array'
%             fn = 'mz';
%           case 'intensity array'
%             fn = 'intensity';
%           case 'zlib compression'
%             compression = 'zlib';
%           case '64-bit float'
%             datatype = 'double';
%           otherwise
%             error('parameter not recognized');
%         end
%       end
%       % Read the raw data
%       alldata = thisbinary.getElementsByTagName('binary');
%       data = alldata.item(0).item(0).getData;
%       % Put it all together
%       tmp.(fn) = xmlbinconvert(data,compression,datatype)';
%     end
%     tmp.n_points = length(tmp.mz);
%     if isempty(scan)
%       scan = tmp;
%     else
%       scan(end+1) = tmp;
%     end
%   end
% end
% 
% function indx = findparams(ref,msparams)
%   flag = strcmp(ref,{msparams.spectrumparamsname});
%   if sum(flag) ~= 1
%     error([ref ' not found'])
%   end
%   indx = find(flag);
% end

function out = xmlbinconvert(data64,compression,datatype)
  % Note: after decoding and, if applicable, decompressing, theoretically
  % one should convert from little-endian to the machine's native
  % representation, but my machine is already little-endian. Laziness wins.
  %http://bryanesmith.com/documents/reading-binary-data-mzml/reading-binary-data-mzml-1.html
  
  % Undo the base64 encoding
  encoder = org.apache.commons.codec.binary.Base64;
  data = encoder.decode(uint8(char(data64)));
  % The remainder was modified from "dunzip",
  % http://www.mathworks.com/matlabcentral/fileexchange/8899-rapid-lossless-data-compression-of-numerical-or-string-variables
  if strcmp(compression,'zlib')
    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
    a=java.io.ByteArrayInputStream(data);
    b=java.util.zip.InflaterInputStream(a);
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
    c = java.io.ByteArrayOutputStream;
    isc.copyStream(b,c);
    out=typecast(c.toByteArray,datatype);
  else
    % untested
    out = typecast(data,datatype);
  end
end
