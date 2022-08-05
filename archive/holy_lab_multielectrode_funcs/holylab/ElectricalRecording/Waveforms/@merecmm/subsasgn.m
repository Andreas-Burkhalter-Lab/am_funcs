function memm = subsasgn(memm,s,b)
% MERECMM/SUBSASGN: assigning values to memmapped MEREC files
% For now this only supports
%    memm.tovolts = false/true
%    memm.blocksize = ...
  
% Copyright 2006 by Timothy E. Holy
  
  switch s.type
   case '.'
    switch s.subs
     case 'blocksize'
       blocksz_in_scans = floor(b/sizeof(memm.type)/ ...
         memm.header.numch);
       if (blocksz_in_scans == 0)
         error('blocksize is too small for any values to be read');
       end
       memm.blocksz_in_scans = blocksz_in_scans;
       % Note the memmap buffer size will be adjusted as needed in subsref
     case 'tovolts'
      memm.tovolts = b;
     case 'contiguous'
      memm.contiguous = b;
     otherwise
      error(['Field ' s.subs{1} ' not recognized']);
    end
  end
  