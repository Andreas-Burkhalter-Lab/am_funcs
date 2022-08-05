function [ljdownindx,ljupindx,hjindx] = whisjclassify(pitch,dpitch)
% WHISJCLASSIFY: find the jump (step) transitions in mouse song syllables
% Syntax:
%   [ljdownindx,ljupindx,hjindx] = whisjclassify(pitch,dpitch)
% where
%   pitch is the dominant frequency vs. time (last point omitted, to give
%     it the same length as dpitch)
%   dpitch is diff(pitch)
% and
%   ljdownindx is the temporal index for all the downward low-jump
%     transitions;
%   ljupindx is the temporal index for all the upward low-jump
%     transitions;
%   hjindx is the temporal index for the high-jump transitions.
% Alternatively, you can call this as
%   hpoly = whisjclassify
% and it will draw polygons corresponding to the cluster definitions.

% Copyright 2005 by Timothy E. Holy
  
  dclust = [65 50; 95 60; 95 40; 65 30; 58 30]*1000;
  %dclust = [65 30; 65 35; 90 46; 90 40]*1000;  % This snips out "line"
  uclust = [45 90; 60 90; 40 55; 30 55; 30 65]*1000;
  %uclust = [30 65; 45 90; 50 90; 30 58]*1000;  % Snips out line
  hclust = [88 77; 88 61; 70 53; 70 64]*1000;

  if (nargin > 0)
    pitchnext = pitch+dpitch;
    ljdown = inpolygon(pitch,pitchnext,dclust(:,1),dclust(:,2));
    ljup = inpolygon(pitch,pitchnext,uclust(:,1),uclust(:,2));
    hj = inpolygon(pitch,pitchnext,hclust(:,1),hclust(:,2));
    ljdownindx = find(ljdown);
    ljupindx = find(ljup);
    hjindx = find(hj);
  else
    % Instead, plot the definitions
    indx = [1:size(uclust,1) 1];
    ljdownindx(1) = line(uclust(indx,1)/1000,uclust(indx,2)/1000,'Color','k');
    indx = [1:size(dclust,1) 1];
    ljdownindx(2) = line(dclust(indx,1)/1000,dclust(indx,2)/1000,'Color','k');
    indx = [1:size(hclust,1) 1];
    ljdownindx(3) = line(hclust(indx,1)/1000,hclust(indx,2)/1000,'Color','k');
  end