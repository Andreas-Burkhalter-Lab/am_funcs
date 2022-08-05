function [out] = fieldout(structure,fieldname,opt);

% function [out] = fieldout(structure,fieldname,type)
%
% routine to extract the information from a particular field
% of a structure array and output it in matrix form (if it needs 
% to be output as a cell array, the optional flag "type = 'cell'" 
% can be used)
%
% See also: subfieldout

if nargin < 3
    opt = struct;
end
if isstr(opt)
    if strcmp(opt,'cell')
        opt = struct;
        opt.cell = 1;
    else
        opt = struct;
        opt.cell = 0;
    end
end
opt = default(opt,'cell',0);
opt = default(opt,'nOut',1);

if ~(length(opt.nOut)==1 & (sum(opt.nOut) == 1))
    sizeOut = opt.nOut;
    out = getfield(structure,char(fieldname));
    return
else
    sizeOut = size(structure);
end
nOut = prod(sizeOut);
if opt.cell
  out = cell(sizeOut);
  for nthEntry = 1:nOut
    out{nthEntry} = getfield(structure,{nthEntry},fieldname);
  end
else
    out = NaNs(sizeOut);
    for nthEntry = 1:nOut
        out(nthEntry) = getfield(structure,{nthEntry},char(fieldname));
    end
end

