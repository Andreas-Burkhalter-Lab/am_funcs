%%% adapted from redo_dFF writtin by Moi Arriaga for Arriaga and Han 2017;
%%% instead of computing dFF after subtracting neuropil, get dff of the neuropil and dff of uncorrected ROI,
%%%  then subtract alpha*neuropil dff from roi dff; 
% using the Arriaga and Han method resulted in negative neuropil-corrected dFF values
%       because mean roi_val minus mean neuropil tended to be negative for any given window length
%       even though roi Ca transients went above neuropil val, so denominator F in dF/F was negative
% new method follows kamani et al; alpha = .7 for kamani, .85 for moi
function dFF_roi_minus_alpha_dFF_neuropil=redo_dFF_AM(F,Fs,varargin) % AM: Fs = imageRate (see dat.ops.imageRate in _proc file)

% alpha = scalefactor for dFF neuropil before subtracting from dFF roi...
%   ...0.7 alpha value taken from Kamani et al. 2016: 'Cooperative subnetworks of molecularly-similar interneurons in mouse neocortex' 

if nargin>2
    time=varargin{1};
else time=300; % AM: ? time = length of window in seconds? Fs is in hz and gets multiplied by 'time' to generate 'window'
end
if nargin>3
    nF=varargin{2}; % AM: nF = neuropil F traces (rows = timepoints, neuropil rois)
end
if nargin>4
    plotOpt =varargin{3};
else plotOpt=0;
end
if nargin>5
    neuropil_alpha_val = 0.7; 
else neuropil_alpha_val = 0.7; 
end

nrois = size(F,2);
wbar = waitbar(0,'Getting dF/F...');
%% roi dF/F
dFF_roi=zeros(size(F));
for iroi=1:nrois % AM: for each ROI (columns of F)...
    F_temp = F(:,iroi);
    numframes=length(F_temp);
    % AM: ? time = length of window in seconds? Fs is in hz and gets multiplied by 'time' to generate 'window'
    % AM: window number of frames are taken before and after each frame
    window=round(Fs*time); 

        F2=zeros(size(F_temp)); %% AM: F2 will contain the 8th F-intensity percentile within the window for each timepoint
        for k=1:length(F_temp) % AM: for each timepoint (rows of nF or of F)...
            % AM: max(1,k-window) starts cut at either frame 1 or window frames before the current timepoint k
            % AM: min(numframes,k+window) takes either the last frame or window frames after the current timepoint k 
            cut=F_temp(max(1,k-window):min(numframes,k+window)); 
            cutsort=sort(cut); %% 
            a=round(length(cut)*.08); % AM: a = index of the 8th percentile of cutsort
            F2(k)=cutsort(a); % AM: for this timepoint, save into F2 the 8th F intensity percentile timepoint within the window
        end
        x=1:numframes;
        if ~isnan(mean(F_temp))
        xBase=F2;
        % AM: get deltaF/F by subtracting 8th percentile of normed trace from normed trace( (nF)-(xBase) == deltaF),...
        %       then divide by 8th percentile of normed trace (xBase == F)
        dFF_roi(:,iroi)=((F_temp)-(xBase))./xBase; 
        
       else
            dFF_roi(:,iroi)=nan(size(F_temp)); % AM: if there are nans, set dFF to nan
        end
    if isvalid(wbar)
        waitbar(iroi/nrois/2,wbar);
    end
end


%% neuropil dF/F
dFF_neuropil=zeros(size(F));
for iroi=1:nrois % AM: for each ROI (columns of F)...
    nF_temp = nF(:,iroi);
    numframes=length(nF_temp);
    % AM: ? time = length of window in seconds? Fs is in hz and gets multiplied by 'time' to generate 'window'
    % AM: window number of frames are taken before and after each frame
    window=round(Fs*time); 

        nF2=zeros(size(nF_temp)); %% AM: nF2 will contain the 8th F-intensity percentile within the window for each timepoint
        for k=1:length(nF_temp) % AM: for each timepoint (rows of nF or of F)...
            % AM: max(1,k-window) starts cut at either frame 1 or window frames before the current timepoint k
            % AM: min(numframes,k+window) takes either the last frame or window frames after the current timepoint k 
            cut=nF_temp(max(1,k-window):min(numframes,k+window)); 
            cutsort=sort(cut); %% 
            a=round(length(cut)*.08); % AM: a = index of the 8th percentile of cutsort
            nF2(k)=cutsort(a); % AM: for this timepoint, save into nF2 the 8th F intensity percentile timepoint within the window
        end
        x=1:numframes;
        if ~isnan(mean(nF_temp))
        xBase=nF2;
        % AM: get deltaF/F by subtracting 8th percentile of normed trace from normed trace( (nF)-(xBase) == deltaF),...
        %       then divide by 8th percentile of normed trace (xBase == F)
        dFF_neuropil(:,iroi)=((nF_temp)-(xBase))./xBase; 
        
       else
            dFF_neuropil(:,iroi)=nan(size(F(:,iroi))); % AM: if there are nans, set dFF to nan
        end

        if plotOpt
            figure
            nSq=ceil(sqrt(size(F,2)));
            for j=1:size(F,2)
            subplot(nSq,nSq,j)
            plot(dFF(:,j));
            end
        end
    if isvalid(wbar)
        waitbar(0.5 + iroi/nrois/2, wbar)
    end

end
if isvalid(wbar)
    close(wbar)
end

%% output
dFF_roi_minus_alpha_dFF_neuropil = dFF_roi - neuropil_alpha_val*dFF_neuropil;