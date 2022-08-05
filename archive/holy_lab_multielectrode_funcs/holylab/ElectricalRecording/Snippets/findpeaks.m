function [tpeak,xmax,samplePeaks,reboundStats] = findpeaks(wave,chanIndex,thresh,left,right,soptions,reboundStats)

% FINDPEAKS: find peaks in wave and return rough times, find times and peak values
% 
% Syntax:
%    [tpeak,xmax,samplePeaks] = findpeaks(wave,chanIndex,thresh,left,right,soptions)
%  
% pre:
%    thresh: a 2x(#ch) matrix, each col is for each ch. row 1 is for
%            polarity==-1 and row 2 is for polarity==1
% 
% post:
%    tpeak: the TIME values on peak points (the prefix t means time). These
%      are unity-offset.
%    xmax: fine times from interpolation, 0 when soptions.interptimes is false;
%    samplePeaks: the SAMPLE values when time=tpeak (not the y values on
%                 interpolated wave)   
%    reboundStats: information on what peaks were thrown out because they
%                  appeared to be the rebound of larger spikes

% History:
%   2004-08-10: 1. Mechanism for screening out 'rebound spikes' added
%               (RCH)
%   ?: modified by Jason for the new snippet header&body format
% 
% Notes:
%    These comments are added by Jason, and they may not be accurate.
% 

nchannels = length(chanIndex);
npts = size(wave,2);
xmax = cell(1,length(chanIndex));
if (soptions.reboundS==0 | soptions.reboundC==0)      % Check to see if going to screen
    reboundFlag = 0;                                  %    for 'rebound' spikes or not
elseif soptions.polarity==-1                           % Only apply rebound screen to
    reboundFlag = 0;                                  %    positive spikes (change later?)
else    
    reboundFlag = 1;
end
for i = 1:nchannels
  w = double(wave(chanIndex(i),:));
  if (soptions.polarity == 1 & reboundFlag==0) % If screening for rebounds, need undershoots
    ihigh = find(w > thresh(2,i));             %    too (will throw out unwanted spikes at end)
  elseif (soptions.polarity == -1)                  
    ihigh = find(w < thresh(1,i));
  else
    % aw = abs(w);
    ihigh = find(w > thresh(2,i) | w < thresh(1,i) );
  end
  % Eliminate values that will cause trouble
  left1 = max([left 1]);
  right1 = max([right 1]);
  ibad = find(ihigh <= left1 | ihigh > length(w)-right1);
  ihigh(ibad) = [];
  % Look for 3-pt extrema
  if (soptions.polarity == 1 & reboundFlag==0)
    whigh = w(ihigh);
    ii3pts = find(whigh > w(ihigh-1) & whigh >= w(ihigh+1));
  elseif (soptions.polarity == -1)
    whigh = w(ihigh);
    ii3pts = find(whigh < w(ihigh-1) & whigh <= w(ihigh+1));
  else
    whigh = w(ihigh);
    ii3pts = find((whigh > w(ihigh-1) & whigh >= w(ihigh+1)) | (whigh < w(ihigh-1) & whigh <= w(ihigh+1)));
  end
  tp1 = ihigh(ii3pts);
  % now tp1 holds peaks' indices (i.e. time or scan#) pointing to items
  %     in w (or aw)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Screen for spikes to toss:
  
  % I. amalgamate those peaks that are
  %  (a) close together
  [tac,indx] = autocorrspike(tp1,soptions.close);
  %  (b) of the same sign
  mnthrsh = mean(thresh(:,i));
  ids = find((w(tp1(indx(1,:)))-mnthrsh) .* (w(tp1(indx(2,:)))-mnthrsh) < 0);
  indx(:,ids) = [];
  %  and (c) do not have a sufficiently-deep trough
  killflag = zeros(size(indx));
  for j = 1:size(indx,2)
    awgap = abs(w(tp1(indx(1,j)):tp1(indx(2,j))) - mnthrsh); % fixed on 9/19/2003
    [mnpk,mnindx] = min(awgap([1 end]));
    if ~any(awgap < soptions.troughdepth*mnpk)
      killflag(mnindx,j) = 1;
    end
  end
  iikill = find(killflag);
  ikill = unique(indx(iikill));
  tp1(ikill) = [];
  
  % II. Toss rebounds of large spikes (must meet all of following criteria
  %   to be tossed):  1) meet criteria for closeness (set by reboundS in 
  %                      snipoptions) 
  %                   2) preceeded by (a) a peak greater than 3x
  %                      its amplitude, and (b) an undershoot greater
  %                      than 2x its amplitude
  %                   3) not followed an undershoot themselves
  %                   4) too low a curvature value (set by reboundC in
  %                      snipoptions)
  if reboundFlag==1
     [tac,indx] = autocorrspike(tp1,soptions.reboundS);   % check (1) closeness 
     triads = [];                                         % check (2) last of triad of high peak & low undershoot
     count=1;
     for in=2:(size(indx,2)-1)
         if (indx(1,in)==indx(1,(in-1)) & indx(2,in)==indx(2,(in+1)))
             triads(:,count)=[indx(1,in); indx(2,(in-1)); indx(2,in)];
             count=count+1;
         end
     end
     if size(triads,2)>0
         notAltPolarity = find( ~( ((w(tp1(triads(3,:)))-mnthrsh).*(w(tp1(triads(2,:)))-mnthrsh)<0) & ...
             ((w(tp1(triads(3,:)))-mnthrsh).*(w(tp1(triads(1,:)))-mnthrsh)>0) ) );
         triads(:,notAltPolarity) = [];
         inohighprev = find( (w(tp1(triads(3,:)))-mnthrsh) > 3.*(w(tp1(triads(1,:)))-mnthrsh) );
         triads(:,inohighprev) = [];
         inolowprev = find( (w(tp1(triads(3,:)))-mnthrsh) > -2.*(w(tp1(triads(2,:)))-mnthrsh) );
         triads(:,inolowprev) = [];  
         for n=1:size(triads,2)                           % check (3) no undershoot
             in = find( indx(1,:)==triads(3,n) );
             if ( (w(tp1(indx(2,in)))-mnthrsh < 0) & (tac(in) < soptions.close) )
                 triads(:,n)=0;
             end
         end
         ihasundershoot = find( triads(1,:)==0 );
         triads(:,ihasundershoot) = [];
     end
     if size(triads,2)>0                                  % check (4) curvature
         v=w(tp1(triads(3,:)))-mnthrsh;
         vp=w(tp1(triads(3,:))+3)-mnthrsh;
         vm=w(tp1(triads(3,:))-3)-mnthrsh;
         curvatures= -( (vp+vm-2.*v) / (2*(3^2)).*v ); 
         prevL = size(reboundStats,2);
         reboundStats(1,(prevL+1):(prevL+size(triads,2))) = ...  % save data to be used as output report
             curvatures;
         reboundStats(2,(prevL+1):(prevL+size(triads,2))) = ...
             tp1(triads(3,:))-tp1(triads(1,:));
         reboundStats(3,(prevL+1):(prevL+size(triads,2))) = ...
             (curvatures > soptions.reboundC);
         reboundStats(4:6,(prevL+1):(prevL+size(triads,2))) = ...
             [tp1(triads(1,:)); tp1(triads(2,:)); tp1(triads(3,:))];
         reboundStats(7,(prevL+1):(prevL+size(triads,2))) = 1;  % flag to show scan numbers haven't yet been adjusted to incluce scanstart
         ihighcurv = find(curvatures > soptions.reboundC);
         triads(:,ihighcurv) = [];
     end
     if size(triads,2)>0
          tp1(triads(3,:)) = [];
     end
  end

  % III. If polarity set to 1 but both signs checked for rebound screen,
  % now eliminate the negative spikes
  if (soptions.polarity==1 & reboundFlag==1)
      ineg = find( w(tp1) < mnthrsh );
      tp1(ineg) = [];
  end
    
  % IV. Eliminate close triggers of opposite sign (if taking samples of both
  %   polarities)  
  if (soptions.polarity==0)
    [tac,indx] = autocorrspike(tp1,soptions.peaktrough);
    iop = find( (w(tp1(indx(1,:)))-mnthrsh) .* (w(tp1(indx(2,:)))-mnthrsh) < 0); % fixed on 9/19/2003
    tp1(indx(2,iop)) = []; % fixed on 9/19/2003 %%TEMP COMMENT OUT
  end
  
  %%%% end screen for spikes to toss
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % if (~isempty(ihigh)) % in version 200303.270
  if (soptions.interptimes & ~isempty(ihigh))
    yminus = w(ihigh(ii3pts)-1);
    yo = w(ihigh(ii3pts));
    yplus = w(ihigh(ii3pts)+1);
    
    xmax{i} = (yminus - yplus)./(2*yplus - 4*yo + 2*yminus);
  else
    xmax{i} = 0;
  end
  tpeak{i} = tp1; 
   % get one more info than old findpeaks() did:
  samplePeaks{i} = w(tp1);
  
end

% Any multichannel stuff could go here

