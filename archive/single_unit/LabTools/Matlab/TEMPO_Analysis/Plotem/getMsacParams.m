function [out] = getMsacParams(msac,pemt,pemh,pemv)

% function [out] = getMsacParams(msac,pemt,pemh,pemv)
% 
% Function to extract parameters (direction and amplitude) for fixational 
% saccades.
%
% INPUTS: msac, pemt, pemh, and pemv.  These are all matrices generated
%         by plotem.m and contain, respectively, the onset and offsets
%         of saccades, the time (from fixation point onset), the horizontal
%         and vertical eye positions.  For more information see plotem.m
%         
% OUTPUTS: The matrix returned by this function contains one row per
%          saccade, and the columns are as follows: 
%          trialnumber, 
%          time of saccade onset,
%          time of saccade offset,
%          direction, 
%          amplitude.
%
% Written by GDLH, 12/19/98.

ntrials = size(msac,1);

% Rotating everything so that columns are trials.
% This makes indexing into the matrices easier.

msac = msac';
pemt = pemt';
pemh = pemh';
pemv = pemv';
startidx = (msac == 1);
endidx = (msac == -1);
starttimes = pemt(startidx);
endtimes = pemt(endidx);
Hcomp = pemh(endidx) - pemh(startidx);
Vcomp = pemv(endidx) - pemv(startidx);
dir = atan2(Vcomp,Hcomp);
amp = sqrt(Hcomp.^2+Vcomp.^2);

nsacs = sum(msac==1);  % Number of saccades in each trial
trialvect = [];
for i = 1:ntrials
  trialvect = [trialvect;i*ones(nsacs(i),1)];
end

out = [trialvect, starttimes, endtimes, dir, amp];
% out = getMsacParams(msac,pemt,pemh,pemv); figure; [u,v] = pol2cart(out(:,3),out(:,4)); compass(u,v); print -dps2; close;