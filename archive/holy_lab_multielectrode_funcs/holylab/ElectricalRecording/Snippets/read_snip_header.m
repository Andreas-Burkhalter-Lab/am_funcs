function header = read_snip_header(file)
% READ_SNIP_HEADER: read snippet header into a matlab structure
%                    
% Syntax:
%    header = read_snip_header(file)
%    
% pre:   
%    file: file name
%    
% post:
%    header: a matlab structure holds all snippet header info
%  
% see also: STRINGIZE_SNIP_HEADER, UPDATE_HEADER 
  
  
   % assume caller is sure the file is snippet file, otherwise here we
   % should test if the file is snippet file by checking magic number at
   % the beginning of the file
  
   header=read_merec_header(file);
   tstrHeader=header.wholeheader;  
   header.snipbeginoffset=str2num(key2value(tstrHeader,'snipbeginoffset'));
   header.snipendoffset  =str2num(key2value(tstrHeader,'snipendoffset'));
   
   tThresholds           =split_dbl(key2value(tstrHeader,'thresh'), ' ');
   header.thresh         =zeros(2, size(tThresholds,2)/2); % make room
                                                           % for row assignment
   header.thresh(1,:)    =tThresholds(1:size(tThresholds,2)/2);
   header.thresh(2,:)    =tThresholds(size(tThresholds,2)/2+1:end); % @todo: error detection
   
   header.numofsnips     =split_int(key2value(tstrHeader,'numofsnips'), ' ');
   header.timesfpos      =split_int(key2value(tstrHeader,'timesfpos'), ' ');
   header.snipsfpos      =split_int(key2value(tstrHeader,'snipsfpos'), ' ');
   header.sniprange      = [header.snipbeginoffset header.snipendoffset];
   
   % now are new fields if need:
   header.inputfile      =key2value(tstrHeader,'snippet input file');
   header.finetimesfpos  =split_int(key2value(tstrHeader,'finetimesfpos'), ' ');
   header.detpeaksfpos   =split_int(key2value(tstrHeader,'detpeaksfpos'), ' ');
   % snippet options:
   header.options=struct('condfilta',   split_dbl(key2value(tstrHeader,'condfilta'), ' '),...
                         'condfiltb',   split_dbl(key2value(tstrHeader,'condfiltb'), ' '),...
                         'detfilt',     makeDetFiltMatrix(key2value(tstrHeader,'detfilt')), ...
                         'polarity',    str2num(key2value(tstrHeader,'polarity')), ...
                         'close',       str2num(key2value(tstrHeader,'close')), ...
                         'troughdepth', str2num(key2value(tstrHeader,'troughdepth')), ...
                         'peaktrough',  str2num(key2value(tstrHeader,'peaktrough')), ...
                         'blocksize',   str2num(key2value(tstrHeader,'blocksize')), ...
                         'interptimes', str2num(key2value(tstrHeader,'interptimes')), ...
                         'interpsnips', str2num(key2value(tstrHeader,'interpsnips'))  ...
                        );
                         

% an aux function:    
function detfiltmatrix=makeDetFiltMatrix(astring)
% pre:
%   astring: ; is used to separated diff column/channel,
%            space is used to separated values within one column/channel,
%            The order of values: the values before first ; are all
%            filter values for the first channels, i.e. filter values are
%            saved column-wise (channel-wise).
   column=split_str(astring, ';');
   detfiltmatrix=[];
   for i=1:size(column, 1) 
      detfiltmatrix(:,i)=split_dbl(column(i,:), ' ');
   end
  
  
