function [fp,fpos] = filehandle(varargin)
% FILEHANDLE: a class for storing information about the disk location of
% a piece of information.
%
% Note this only tells you where to find the file and basic info like the
% byteorder; there must be additional data to specify position in the
% file, and functions to specify how the information is formatted.
%
% Syntax:
%   fh = filehandle('propertyname1',value1,'propertyname2',value2,...)
% where valid propertynames are:
%   abspathstr   (absolute path string, must be terminated by filesep)
%   filename
%   fid
%   machfmt (see FOPEN)
%   uselfs (use Large File System functions, true or false)
%
% It's also valid to say
%   fh = filehandle(s)
% where s is a structure containing fields with the property names.
%
% Alternatively, you can say
%   [fh,fpos] = filehandle(fid)
% to set up the filehandle and get the file position, given a file
% identifier such as that opened by FOPEN.
%
% See also: FOPEN.

% Copyright 2005 by Timothy E. Holy

  if (nargin == 1 && isa(varargin{1},'filehandle'))
    % Copy constructor
    fp = varargin{1};
    return;
  end
  
  fp.abspathstr = '';
  fp.filename = '';
  fp.fid = [];
  fp.machfmt = 'n';
  fp.uselfs = 0;

  if (nargin == 1)
    if isnumeric(varargin{1})
      % Initialize from a file identifier
      [fp.filename,mode,fp.machfmt] = fopen(varargin{1});
      % Map machine format to standard single-character names
      % (this facilitates using the "isequal" command in comparing real &
      % faked saves with imwrite & moviewrite)
      longnames = {'ieee-le','ieee-be','ieee-le.l64','ieee-be.l64'};
      shortnames = {'l','b','a','s'};
      indx = strmatch(fp.machfmt,longnames,'exact');
      if ~isempty(indx)
        fp.machfmt = shortnames{indx};
      end
      fpos = ftell(varargin{1});
    elseif isstruct(varargin{1})
      % Initialize from a structure
      fp = copyfields(varargin{1},fieldnames(fp),fp);
    else
      error('Constructor syntax not recognized');
    end
  else
    propnames = varargin(1:2:end);
    values = varargin(2:2:end);
    extranames = setdiff(propnames,fieldnames(fp));
    if ~isempty(extranames)
      errmsg = [sprintf('Unknown propertynames:\n') ...
                sprintf('%s\n',extranames{:})];
      error(errmsg)
    end
    for i = 1:length(propnames)
      fp.(propnames{i}) = values{i};
    end
  end
  
  fp = class(fp,'filehandle');
