%%%%% AM comments 18/5/17

function dFF=redo_dFF(F,Fs,varargin) % AM: Fs = imageRate (see dat.ops.imageRate in _proc file)
if nargin>2
    time=varargin{1};
else time=300; % AM: ? time = length of window in seconds? Fs is in hz and gets multiplied by 'time' to generate 'window'
end
if nargin>3
    nF=varargin{2}; % AM: nF = neuropil F traces (rows = timepoints, neuropil rois)
end
if nargin>4
    plotOpt =1;
    figure
else plotOpt=0;
end
    
% F=mouse(3).F0{14}.raw{4};

dFF=zeros(size(F));
for j=1:size(F,2) % AM: for each ROI (columns of F)...
junk=F(:,j)-.85*nF(:,j); % AM: junk = for each timepoint for this roi, subtract 85% of corresponding neuropil trace at that timepoint
numframes=length(junk);
% AM: ? time = length of window in seconds? Fs is in hz and gets multiplied by 'time' to generate 'window'
% AM: window number of frames are taken before and after each frame
window=round(Fs*time); 

        junk2=zeros(size(junk)); %% AM: junk2 will contain the 8th F-intensity percentile within the window for each timepoint
        for k=1:length(junk) % AM: for each timepoint (rows of junk or of F)...
            % AM: max(1,k-window) starts cut at either frame 1 or window frames before the current timepoint k
            % AM: min(numframes,k+window) takes either the last frame or window frames after the current timepoint k 
            cut=junk(max(1,k-window):min(numframes,k+window)); 
            cutsort=sort(cut); %% AM: cutsort = sort frames in the window by pixel pixel intensity (normed to F(:,j)-.85*nF(:,j))
            a=round(length(cut)*.08); % AM: a = index of the 8th percentile of cutsort
            junk2(k)=cutsort(a); % AM: for this timepoint, save into junk2 the 8th F intensity percentile timepoint within the window
        end
        x=1:numframes;
        if ~isnan(mean(junk))
%         fitData=fit((x)',junk2,'exp2'); % AM: for this ROI, fit 2-term exponential model to the 8th percentile neuropil-normalized trace in junk2 (AM commented out b/c not used)
%         xBase=((fitData.a)*exp((fitData.b)*(x))+(fitData.c)*exp((fitData.d)*(x)))';
        xBase=junk2;
        % AM: get deltaF/F by subtracting 8th percentile of normed trace from normed trace( (junk)-(xBase) == deltaF),...
        %       then divide by 8th percentile of normed trace (xBase == F)
        dFF(:,j)=((junk)-(xBase))./xBase; 
        
               
        
        else
            dFF(:,j)=nan(size(F(:,j))); % AM: if there are nans, set dFF to nan
        end
        if plotOpt
            nSq=ceil(sqrt(size(F,2)));
            subplot(nSq,nSq,j)
            plot(junk);
            hold on;
            plot(xBase);
        end
end
        if plotOpt
            figure
            nSq=ceil(sqrt(size(F,2)));
            for j=1:size(F,2)
            subplot(nSq,nSq,j)
            plot(dFF(:,j));
            end
        end

% Old dFF method
% if nargin>2
%     time=varargin{1};f
% else time=15;
% end
% 
% dFF=zeros(size(F));
% for j=1:size(F,2)
% junk=F(:,j);
% %         junk=junk;
% numframes=length(junk);
% % numwindow=numframes/(Fs*15);
% % % numwindows=50;
% %  window=round(numframes/numwindow);
% window=round(Fs*time);
%         junk2=zeros(size(junk));
%         for k=1:length(junk)
%             cut=junk(max(1,k-window):min(numframes,k+window));
%             cutsort=sort(cut);
%             a=round(length(cut)*.08);
%             junk2(k)=cutsort(a);
%         end
%         dFF(:,j)=(junk./junk2);
%         maxval=max(dFF(:,j));
%         dFF(:,j)=(dFF(:,j)-median(dFF(:,j)))/max((dFF(:,j)-median(dFF(:,j))));
%         dFF(:,j)=maxval*dFF(:,j);
% end