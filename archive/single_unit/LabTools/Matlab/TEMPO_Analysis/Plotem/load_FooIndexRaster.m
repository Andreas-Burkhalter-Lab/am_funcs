function [foo, paradigm, index, raster] = load_FooIndexRaster(efile,path)
% [foo, index, raster, paradigm] = load_FooIndexRaster(efile,path)
% This function loads the .s, .ras, and .index files that are created by the
% spxxx "spikes" program.  I have done a major overhaul of this now so that it
% is no longer dependent on various Unix commands.  It's much simpler now.

% Note that if the function is called with only 1 (or 2) output argument(s), then only
% the .s (and .index) files will be loaded

% I added code to return the paradigm # in the last output argument
% Greg DeAngelis, 5/11/99

sfile = [path '\\' efile '.snew'];
if fopen(sfile)==-1
   sfile = [path '\\' efile '.s'];
end
[header, foo] = hdrload(sfile);


% here, I extract the paradigm # from the header of the .s file
line = header(3,:);
whitespace = find(line == 32);
paradigm = str2num(line(whitespace(1)+1:whitespace(2)-1));
%foo

if (nargout > 2)	%if we want the index file too
   indexfile = [path '\\' efile '.index'];
   [header, index] = hdrload(indexfile);
   % append a column of zeros to the end of index.  This column is the offse
   % for FixPt On, since that's how the values are set in spikes. 
	% (MNS 10/11/94)  Kept in by GCD 5/11/99
   index = [index zeros(size(index,1),1)]; 
   %index
end
   
if (nargout > 3)	%if we want the rasters
   rasterfile = [path '\\' efile '.ras'];
   [header, raster] = hdrload(rasterfile);
   %raster
end

