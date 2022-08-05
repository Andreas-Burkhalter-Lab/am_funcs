function data = cellsort_by_celltimes(data,inname,outname,celldimension);

% Will take the information in the 'inname' field and, assuming its
% indexing matches that of stiptimes, break it up based on celltimes'
% sorting decisions.  Very inefficient - just blunt force.

if nargin < 4
    celldimension = 1;
end
if celldimension ~= 1
    error('Oops, haven''t implimented this yet.  Sorry.  Prob need to use repmat?')
end

nIntervals = length(data);
for nthInterval = 1:nIntervals
    outinfo = [];
    dat = data(nthInterval);
    ininfo = getfield(data,{nthInterval},inname);
    if (iscell(ininfo) & length(ininfo)==1)
        ininfo = ininfo{1};
    end
    nCells = length(dat.cellnums);
    for nthCell = 1:nCells
        celltimes = dat.celltimes{nthCell};
        sniptimes = dat.sniptimes{1};
        i_cell = findainb(celltimes,sniptimes);
        outinfo{nthCell} = ininfo(i_cell,:);
    end
    data = setfield(data,{nthInterval},outname,outinfo);
end

