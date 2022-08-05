function ipout = imphyscopy(ipin,fields,ipdest)
% IMPHYSCOPY: copy fields between imphys structures
%
% Syntax:
%   ipout = imphyscopy(ipsource,fields)
%   ipout = imphyscopy(ipsource,fields,'check')
%   ipout = imphyscopy(ipsource,fields,ipdest)
% where
%   ipsource, ipdest, and ipout are imphys/imphysp structure arrays;
%   fields is a set of field names. You may also use the following
%     category abbreviations:
%       'all': all fields
%       'public': all but 'private' fields (see IMPHYS)
%       'dimension', 'intensity', 'storage', 'timing', 'info',
%         'analysis', 'private': as described in IMPHYS and IMPHYSP.
%
% The first syntax is used when you want to create a new imphys structure
% from a subset of fields in another imphys structure.
%  
% The second syntax is used when you have an array of imphys structures,
% and you want to assert that the values of the stated fields are
% identical across different structure elements.  Equality is checked,
% and an error is generated if it proves to be false.  On output, ipout
% contains a single element containing the (consensus) field values.
%
%  The third syntax is used when you want to copy data from one imphys
%  structure (ipsource) over to another imphys structure (ipdest).  The
%  precise behavior depends on the sizes of ipsource and ipdest:
%    1. If ipsource and ipdest are both scalar structures, or are
%       structure arrays with the same number of elements, the copy is
%       performed directly.
%    2. If ipsource is an array of structures and ipdest is a scalar,
%       then it is first checked that all the elements of ipsource are
%       identical in the chosen fields; if not, an error is issued.
%
% See also: IMPHYS, IMPHYSP.
  
% Copyright 2005 by Timothy E. Holy
  
  if ~iscell(fields)
    fields = {fields};
  end
  nfields = length(fields);
  expfields = {};
  catfields = {'dimension',{'xrange','x2um','yrange',...
                            'y2um','zrange','z2um','stackz'};...
               'intensity',{'imrange'}; ...
               'timing',{'stacknum','stacktime'}; ...
               'info',{'tag','date','stimulus'}; ...
               'storage',{'imfile','headerfile','stackfpos','imfilefmt',...
                          'width','height','depth','pixelorder',...
                          'imfilemachfmt','imfileprec','tform'}; ...
               'analysis',{'oimfile','ostacknum','ostackweight', ...
                           'filterwidth','filtersize','resampmag', ...
                           'resampmethod'}; ...
               'private',{'image'}};
  ncategories = size(catfields,1);
  privateindx = strmatch('private',catfields(:,1),'exact');
  publicindx = setdiff(1:ncategories,privateindx);
  catfields(end+1,1:2) = {'public',cat(2,catfields{publicindx,2})};
  catfields(end+1,1:2) = {'all',cat(2,catfields{1:ncategories,2})};
  for i = 1:nfields
    indx = strmatch(fields{i},catfields(:,1),'exact');
    if ~isempty(indx)
      expfields{end+1} = catfields{indx,2};
    else
      expfields{end+1} = fields{i};
    end
  end
  expfields = cat(2,expfields{:});
  % Take care of 2nd syntax
  mode = 'copy';
  if (nargin > 2)
    if isstruct(ipdest)
      ipout = ipdest;
      if (length(ipin) > length(ipdest))
        if (length(ipdest) > 1)
          error('ipdest is of the wrong size');
        end
        mode = 'checkcopy';
      end
    elseif (isstr(ipdest) & strcmp(ipdest,'check'))
      mode = 'check';
      if (nargout > 0)
        mode = 'checkcopy';
      end
    else
      error('Syntax not recognized');
    end
  end
  if strcmp(mode,'copy')
    % Do the copying
    for i = 1:length(ipin)
      for j = 1:length(expfields)
        if isfield(ipin,expfields{j})
          ipout(i).(expfields{j}) = ipin(i).(expfields{j});
        end
      end
    end
  else
    % Check, and possibly copy
    for j = 1:length(expfields)
      if isfield(ipin,expfields{j})
        tmp = {ipin.(expfields{j})};
        if ~isequal(tmp{:})
          error('Not all fields are filled with the same values');
        end
        if strcmp(mode,'checkcopy')
          ipout.(expfields{j}) = ipin(1).(expfields{j});
        end
      end
    end
  end

  