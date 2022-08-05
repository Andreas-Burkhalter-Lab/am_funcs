function sh = snipfile2sortheader(filenames)
% SNIPFILE2SORTHEADER: create a sortheader from a snippet file
%
% This function forms a gateway between the file format used to store
% spike snippets and the autosorting program.  It parses the disk file
% and organizes the relevant information into the fields of a structure.
%
% Syntax:
%   sh = snipfile2sortheader(filenames)
% where
%   filenames is a cell array, each entry containing the name of a
%     snippet file;
% and
%   sh is a structure array of sortheaders.
%
% sortheaders must contain the following fields, used by AUTOSNIP and
% CASS:
%   (To be written---just try it and see!)
%
% See also: AUTOSORT, CASS.
  
% Copyright 2005 by Timothy E. Holy
  
  if ischar(filenames)
    filenames = {filenames};
  end
  nfiles = length(filenames);
  for i = 1:nfiles
    % Determine whether this file is bigger than 2GB (use Large File System)
    uselfs = should_use_lfs(filenames{i});
    % Read in the header, in a way that permits legacy files to be read.
    h = readheader(filenames{i},struct('headertype','snip'));
    % Massage into prettier/more rational format
    machfmt = h.endian;
    htmp = rmfield(h,{'headersize','endian','wholeheader',...
                       'snipbeginoffset','snipendoffset',...
                       'inputfile'});
    [pathname,basename,ext] = fileparts(filenames{i});
    if (isempty(pathname) || ~strmatch(pathname(1),filesep))
      pathname = [pwd filesep pathname];
    end
    htmp.fh = filehandle('abspathstr',pathname,...
                         'filename',[basename ext],...
                         'machfmt',machfmt,...
                         'uselfs',uselfs);
    htmp.polarity = h.options.polarity;
    if (htmp.polarity ~= 0 && length(htmp.thresh) > 1)
      polarity = [-1 1];
      htmp.thresh = htmp.thresh(find(htmp.polarity == polarity),:);
    end
    htmp.timesprec = 'int32';
    htmp.snipsprec = 'int16';
    htmp.detpeaksprec = 'int16';
    sh(i) = htmp;
  end
  % Finally, check to make sure that these files have a unique start time
  % (e.g., if the user accidently supplies a .ssnp and a _clean.ssnp
  % file)
  starttime = sortheader_absolute_starttime(sh);
  if (length(unique(starttime)) < length(starttime))
    error(['Files do not have unique start times---perhaps there has been ' ...
           'duplication?']);
  end