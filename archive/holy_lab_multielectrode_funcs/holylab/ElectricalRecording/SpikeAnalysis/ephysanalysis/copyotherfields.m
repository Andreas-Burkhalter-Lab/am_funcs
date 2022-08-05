function sout = copyotherfields(sin,fieldstoskip,options)
% COPYOTHERFIELDS: copy some fields from one structure to a new one
% sout = copyotherfields(sin,fieldstoskip)
% where
%   sin is a structure array
%   fieldstoskip is a cell array of field names that are to be omitted
%
  fields = fieldnames(sin);
  nfields = length(fields);
  nstructs = prod(size(sin));
  temp = cell(1,nfields);
  for i = 1:nstructs
    for j = 1:nfields
      if strmatch(fields{j},fieldstoskip,'exact')
        temp{j} = {[]};
      else
        temp{j} = {sin(i).(fields{j}) };
      end
    end
    m = [fields';temp];
    sout(i) = struct(m{:});
  end
  sout = reshape(sout,size(sin));

% Here's the old version: it didn't preserve the overall
% form of the structure
%   fields = fieldnames(sin);
%   fieldstocopy = setdiff(fields,fieldstoskip);
%   nfields = length(fieldstocopy);
%   nstructs = prod(size(sin));
%   temp = cell(1,nfields);
%   for i = 1:nstructs
%     for j = 1:nfields
%       temp{j} = {getfield(sin(i),fieldstocopy{j})};
%     end
%     m = [fieldstocopy';temp];
%     sout(i) = struct(m{:});
%   end
%   sout = reshape(sout,size(sin));
