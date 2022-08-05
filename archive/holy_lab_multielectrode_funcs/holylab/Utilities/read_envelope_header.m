function header = read_envelope_header(file)
% READ_ENVELOPE_HEADER: read envelope header into a matlab structure
%                    
% Syntax:
%    header = read_envelope_header(file)
%    
% pre:   
%    file: file name
%    
% post:
%    header: a matlab structure holds all envelope header info
%  
  
   % assume caller is sure the file is envelope file, otherwise here we
   % should test if the file is envelope file by checking magic number at
   % the beginning of the file
  
   header=read_merec_header(file);
   tstrHeader=header.wholeheader;  
   header.decimate=str2num(key2value(tstrHeader,'block size'));
   
   % now are new fields if need:
   header.inputfile      =key2value(tstrHeader,'envelope input file');
   
