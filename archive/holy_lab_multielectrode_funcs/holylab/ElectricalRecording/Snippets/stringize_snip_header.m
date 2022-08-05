function strHeader = stringize_snip_header(header, soptions, srcfilename)
% STRINGIZE_SNIP_HEADER: convert header in matlab structure format to
%                        string and also add snip options to the string
% strHeader = stringize_snip_header(header, soptions)
% 
% The outputs are file position values for items that are likely
%        to need updating later:
%        nsnppos: file positions for the # of snips/channel vector
%        timepospos: file positions for the spiketime-file-positions vector
%        snippospos: file positions for the spikewaveform-file-positions
%                    vector
%                    
%  
%  Notes: 1. This func is used by snippetfile, SubSnipFile,
%         and WriteSnipFile. It may break the latter two.
%         2. This func doesn't write header to file, and it just creates a
%         string header. Please use update_header to write the header
%         content to file.
%         3. This func assumes header.wholeheader field is consistent
%         with other fields, especially headersize field
%         4. todo: need to auto adjust headersize field
%         
%  todo: test this w/ SubSnipFile and WriteSnipFile.
%  
%  limitation: header must be a new ascii format header
%  
%  obsolete: WriteSnipHeader
%  
%  see also: UPDATE_HEADER, STRINGIZE_SNIP_HEADER
  
  
   strHeader=header.wholeheader;
   strHeader=update_magic(strHeader, 'SNIPPET');
   % the string concatenation maybe is inefficient, but it is ok and make
   % code more readable:
   strHeader=[strHeader '[snippet]' char(10)];
   strHeader=[strHeader 'snipbeginoffset=' num2str(header.snipbeginoffset) char(10)]; % use int2str() ?
   strHeader=[strHeader 'snipendoffset=' num2str(header.snipendoffset) char(10)];
   strHeader=[strHeader 'thresh=' num2str(header.thresh(1,:)) ' ' num2str(header.thresh(2,:))  char(10)];
        % @notes: thresholds are saved row-wise
   
   % following 5 fields may be empty:
   numch = length(header.channels); % 
   if ~isfield(header,'numofsnips')
      header.numofsnips = zeros(1, numch);
   end
   if ~isfield(header,'timesfpos')
      header.timesfpos  = zeros(1, numch);
      header.snipsfpos  = zeros(1, numch);
      header.finetimesfpos  = zeros(1, numch);
      header.detpeaksfpos   = zeros(1, numch);
   end
   strHeader=[strHeader 'numofsnips=' int2str(header.numofsnips) char(10)];
   strHeader=[strHeader 'timesfpos=' int2str(header.timesfpos) char(10)];
   strHeader=[strHeader 'snipsfpos=' int2str(header.snipsfpos) char(10)];

   % now are the new fields:
   strHeader=[strHeader 'snippet input file=' srcfilename char(10)];
   strHeader=[strHeader 'finetimesfpos=' int2str(header.finetimesfpos) char(10)];
   strHeader=[strHeader 'detpeaksfpos=' int2str(header.detpeaksfpos) char(10)];
   
   % snippet options:
   strHeader=[strHeader 'condfilta=' int2str(soptions.condfilta) char(10)];
   strHeader=[strHeader 'condfiltb=' int2str(soptions.condfiltb) char(10)];
   strHeader=[strHeader 'detfilt=' stringize_detection_filter(soptions.detfilt) char(10)];
   strHeader=[strHeader 'polarity=' num2str(soptions.polarity) char(10)];
   strHeader=[strHeader 'close=' num2str(soptions.close) char(10)];
   strHeader=[strHeader 'troughdepth=' num2str(soptions.troughdepth) char(10)];
   strHeader=[strHeader 'peaktrough=' num2str(soptions.peaktrough) char(10)];
   strHeader=[strHeader 'blocksize=' num2str(soptions.blocksize) char(10)];
   strHeader=[strHeader 'interptimes=' num2str(soptions.interptimes) char(10)];
   strHeader=[strHeader 'interpsnips=' num2str(soptions.interpsnips) char(10)];
   
   
   return
   
   
   
function strDetFilter=stringize_detection_filter(detfilt)
  
   % 
   if(isempty(detfilt))
      strDetFilter = '';  
   else
      strDetFilter = num2str(detfilt(:,1)');
      for i=2:size(detfilt, 2) 
         strDetFilter=[strDetFilter ';' num2str(detfilt(:,i)')];
      end
   end
   
  
   return 
