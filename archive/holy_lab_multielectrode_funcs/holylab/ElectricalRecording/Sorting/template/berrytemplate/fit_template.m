function fit_template(templateFile, fileToFit, optionsForFit)
% fit_template(templateFile, fileToFit, optionsForFit)
% templateFile: a *.fine_cluster file  
   
   if(nargin==2)
      optionsForFit=struct;
   end
   optionsForFit=fillOptions(optionsForFit);

   templateVars=load(templateFile, '-mat');
   templates=templateVars.templates;
   channels=templateVars.channels;
   thresh=templateVars.thresh;
   medv=templateVars.medv;
   
   rawClusterFile=templateVars.rawClusterFilename;
   tt=load(rawClusterFile, 'options', '-mat');
   rawClusterOptions=tt.options;
   % TODO: raw cluster file may be too big to load
   
   % make the variables required by MultiChannelFitV4(): T and ShiftMatrix
   % NOTE: MultiChannelFitV4() requires its args single precision
   % NOTE: MultiChannelFitV4() assumes 64 samples per snippet
   %  MultiChannelFitV4(t,T,ShiftMatrix,Lambda,NumberOfOverlap,I);
   %    t: waveform; I: 
   nTemplates=length(templates);
   shifts=make_vector(optionsForFit.shifts, 'col');
   nShifts=length(shifts);
   % make ShiftMatrix
   ShiftMatrix=zeros(nTemplates*nShifts+1, 2, 'single'); % NOTE: col 3 and 4 will be all zeros; last row is also all zeros
   for templateIdx=1:nTemplates
      ShiftMatrix([1:nShifts]+(templateIdx-1)*nShifts, 1:2)=[ones(nShifts,1)*templateIdx shifts];
   end
   % make T (i.e. templates with shift)
   tT=cell(1, nTemplates);
   for templateIdx=1:nTemplates
      curTemplate=reshape(templates{templateIdx}, [], length(channels));
      tT{templateIdx}=make_shifts(curTemplate', shifts);
   end % for, each template
   T=cat(3, tT{:});
   T(:,:,end+1)=0;
   
%    T=zeros(size(ShiftMatrix,1), length(templates{1}));
%    snipLen=length(templates{1})/length(channels);
%    padding=zeros(snipLen, length(channels)); % NOTE: b/c data to fit and templates are all normalized, use 0 to pad
%    for templateIdx=1:nTemplates
%       templateWithoutShift=reshape(templates{templateIdx}, snipLen, length(channels));
%       rowFrom=(templateIdx-1)*nShifts;
%       for idxShift=1:nShifts
% 	 curShift=shifts(idxShift);
% 	 if(curShift>=0) 
% 	    % shift to right
% 	    tt=[padding(1:curShift,:); templateWithoutShift(1:end-curShift,:)];
% 	    T(rowFrom+idxShift,:)=tt(:);
% 	 else 
% 	    % shift to left
% 	    curShift=-curShift;
% 	    tt=[templateWithoutShift(curShift+1:end,:); padding(1:curShift,:)];
% 	    T(rowFrom+idxShift,:)=tt(:);
% 	 end
%       end % for, each shift
%    end % for, each template

   % now we can fit data
   for templateIdx=1:nTemplates
      % fitting(templateIdx)=struct('eventTime', [], 'shift', [], 'shiftedTime', []); % in this way, each template
      %    is one element of struct array fitting; every field of a struct is a vector
      fitting{templateIdx}=struct('shiftedTime', {}, 'shift', {}, 'eventTime', {}, 'rowIdxToShiftMatrix', {}); 
      %     in this way, each template (biology cell) is one element of cell array fitting; each matlab cell is a
      %     struct array; each struct represents one spike.
   end % for, each template
   
   % next var is struct array, each struct represents >=1 spike(s) at one time;
   % each field of a struct is a vector (or scalar when only one spike);
   % that is, it is another view of fitting;
   % initialized to 0s, to make code logic easier (will be removed later);
   % Also field idxBlock will be removed later.
   spikeCollection=struct('shiftedTime', 0, 'shift', 0, 'eventTime', 0, 'rowIdxToShiftMatrix', 0, 'templateIdx', 0, 'idxBlock', 0);

   % file to save result
   if(isfield(optionsForFit, 'fileToSave'))
      fileToSave=optionsForFit.fileToSave;
   else
      fileToSave=replace_extension(fileToFit, '.fit');
   end
   
   isSaveRes=optionsForFit.isSaveRes;
   if(isSaveRes)
      resFile=replace_extension(fileToSave, '.residual');
   end
   
   memm = merecmm(fileToFit,'tovolts',true,'contiguous',true);
   header=memm.header;
   blockSize=optionsForFit.blockDuration*header.scanrate; % include overlapping
   overlapSize=optionsForFit.overlapDuration*header.scanrate;
   pureBlockSize=blockSize-overlapSize; % the block size w/o overlapping
   nBlocks=ceil(header.nscans/pureBlockSize);
   for idxBlock=1:nBlocks
      scanNumFrom=(idxBlock-1)*pureBlockSize+1;
      scanNumTo  =min(scanNumFrom+blockSize-1, header.nscans);
      wave=memm(channels, [scanNumFrom scanNumTo]);
      [d,N] = size(wave); % d: #channels; N: #scans
      thresh_mtrx = repmat(thresh,1,N);
      
      allWave=wave-repmat(medv,1,N); % allWave: normalized wave for all times within the block

      % Find the points in time that exceeded the threshold
      absNormalizedWave=abs(allWave);
      isbig = absNormalizedWave > thresh_mtrx; % just big enough
      
      % 3-pts maxima
      if(optionsForFit.isUsePeakCriteria)
         isPeak= absNormalizedWave(:,[1 1:end-1]) < absNormalizedWave  &  absNormalizedWave >= absNormalizedWave(:,[2:end end]);
         isbig=isbig & isPeak; % big peak
      end
      
      if(d>1)
         keep_timeslice = max(isbig); % those times that something happened on at least one electrode
         % TODO: use min(isbig) if requiring "sth happened on all electrodes"
      else
         keep_timeslice = isbig;
      end
      eventTime = find(keep_timeslice); % record the times. NOTE: eventTime is the col indices to allWave
      
      % isMergeAdjacent=0;
      if(optionsForFit.isMergeAdjacent)
	 ibreak = find(diff(eventTime) >= 20); % TODO: similar to cluster_events_raw(). Make it a param.
	 eventTime=eventTime([1 ibreak+1]);
      end
      
      boundaryDistance=32*3; % TODO: make it arg. may need change cluster_events_raw() too.
      % i.e. snippet width * 1.5. The reason is: shift < width/2,
      %      and sniprange(2)<width
      eventTime=eventTime(eventTime>boundaryDistance & eventTime<N-boundaryDistance);
      
      if(isfield(optionsForFit, 'lambda'))
         lambda=optionsForFit.lambda; % NOTE: this option is for debug purpose
      else
         lambda=mean(thresh);   % NOTE: only these two args (lambda and maxNumberOfOverlap) can be double
      end
      

      maxNumberOfOverlap=optionsForFit.maxNumberOfOverlap;

      
      %ST=MultiChannelFitV4(t,             T,        ShiftMatrix,        Lambda,NumberOfOverlap,    I);
      % singleT=single(T); single_eventTime=single(eventTime);
      [ST, res]=MultiChannelFit(allWave,T,ShiftMatrix,lambda,maxNumberOfOverlap, eventTime, ...
         -rawClusterOptions.sniprange(1)+1);
      
      ST=ST(:, [2 1]); % swap the columns
      
      ST=ST(find(ST(:,1)~=0),:); % remove zero-template

      % following checking is not necessary
%       % here assume col 1 of ST is sorted so that for any neoron, fitting{neoron}.eventTime (and .shiftedTime) are added increasingly
%       % TODO: sort it otherwise
%       if(any(sort(ST(:,1))~=ST(:,1))) % OR: if(~isequal(sort(ST(:,1)), ST(:,1)))
% 	 error('not implemented yet');
%       end
      
      % save residual
      if(isSaveRes)
        clear tt
        tt.(['chunk' num2str(idxBlock)])=res;
        % fprintf('now chunk %d\n', idxBlock);
        if(idxBlock>1)
           save(resFile, '-struct', 'tt', '-mat', '-v7.3', '-append');
        else
           save(resFile, '-struct', 'tt', '-mat', '-v7.3');
        end
      end


      for idxMatch=1:size(ST,1)
	 curMatch=ST(idxMatch,:);
	 templateIdx=ShiftMatrix(curMatch(2), 1);
         clear curSpike;
	 curSpike.shiftedTime=curMatch(1)+scanNumFrom-1; % now scan # is relative to the whole file
	 curSpike.shift      =double(ShiftMatrix(curMatch(2), 2));
	 curSpike.eventTime  =curSpike.shiftedTime-curSpike.shift;% shiftedTime is the sum of shift and eventTime
         curSpike.rowIdxToShiftMatrix=curMatch(2); % also dim-3 index to T
	 
	 %fitting(templateIdx).eventTime(end+1)  =curSpike.eventTime;
	 %fitting(templateIdx).shift(end+1)      =curSpike.shift;
	 %fitting(templateIdx).shiftedTime(end+1)=curSpike.shiftedTime;
	 fitting{templateIdx}(end+1)=curSpike;
         
         % now update spikeCollection
         curSpike.templateIdx=templateIdx;
         curSpike.idxBlock=idxBlock;
         if(spikeCollection(end).eventTime(end)==curSpike.eventTime)
            % avoid double counting
            if(spikeCollection(end).idxBlock(end)==curSpike.idxBlock)
               spikeCollection(end)=merge_struct(spikeCollection(end), curSpike);
            % else, not in same block, double counting
            end
         elseif(spikeCollection(end).eventTime(end)<curSpike.eventTime)
            spikeCollection(end+1)=curSpike;
         % else curSpike is double counted, just ignore   
         end
         
      end % for, each row of ST
   end % for, each block
   
   % b/c we use shiftedTime to clean double-count and refractory violation,
   % here we should sort shiftedTime (even if eventTime is in order, +shift may cause shiftedTime out of order)
   for templateIdx=1:nTemplates
      [tt, ttIndices]=sort([fitting{templateIdx}.shiftedTime]);
      fitting{templateIdx}=fitting{templateIdx}(ttIndices);
   end % for, each template
   % get rid of the double counted events due to block-by-block processing
   %  and CleaningRefractoryViolation
   for templateIdx=1:nTemplates
      ibreak=find(diff([fitting{templateIdx}.shiftedTime]) > 1 ); % NOTE: hardcoded 1. Zero is enough for getting rid
                                                                  % of doublecount only.
      if(~isempty(fitting{templateIdx}))                                                           
         fitting{templateIdx}=fitting{templateIdx}([1 ibreak+1]); % +1: say y=diff(x), indiceToY+1 is indicesToX.
      end % if, find at least one fit
   end
   
   % now clean up spikeCollection
   spikeCollection(1)=[]; % get rid of the first one
   % get rid of field idxBlock
   spikeCollection=rmfield(spikeCollection, 'idxBlock');
   
   
   save(fileToSave, 'fileToFit', ...
	'fitting', ...
        'spikeCollection', ...
        'T', ...
        'ShiftMatrix', ...
	'templateFile', ...
	'optionsForFit', ...
	'-mat');

     
function result=merge_struct(s1, s2)
   names=fieldnames(s1);
   for idx=1:length(names)
      name=names{idx};
      result.(name)=[s1.(name) s2.(name)];
   end

     
function options=fillOptions(options)
   if(~isfield(options, 'blockDuration'))
      options.blockDuration=10; % TODO: default to 10 sec
   end
   if(~isfield(options, 'overlapDuration'))
      options.overlapDuration=0.05; % TODO: default to 50ms
   end

   default_options('isUsePeakCriteria', 1);
   % default_options('shifts', -17:48);
   default_options('shifts', -10:10);
   % default_options('maxNumberOfOverlap', 5);
   default_options('maxNumberOfOverlap', 2);
   default_options('isSaveRes', 0);
   
   if(~isfield(options, 'isMergeAdjacent'))
      options.isMergeAdjacent=0; % default to false
   end
   
