function [dist,path,ss,pitchx] = whis_fittemplate(pitch,template,x)
% WHIS_FITTEMPLATE: fit pitch waveforms to a template using DTW
% Syntax:
%   [dist,path,shiftscale] = whis_fittemplate(pitch,template)
%   [dist,path,shiftscale,pitchx] = whis_fittemplate(pitch,template,x)
% where
%   pitch is a cell array of pitch trajectories;
%   template is the waveform to be aligned against;
%   x (optional) contains the x-coordinates of template;
% and
%   dist is a vector giving the sum-squared distance between the warped
%     pitch trajectories and the template;
%   path is a cell array of n-by-2 matrices, giving the DTW path;
%   shiftscale is a nwaveforms-by-2 matrix, where the best-fitting
%     template to pitch{i} is
%         shiftscale(i,2)*template + shiftscale(i,1)
%   pitchx is a cell array of vectors, containing the template's
%     corresponding x-coordinate as a function of time throughout the
%     pitch trajectory.
%
% Copyright 2005 by Timothy E. Holy
  
  nw = length(pitch);
  
  for i = 1:nw
    progress_bar(struct('progress',i,'max',nw));
    cw = pitch{i}; % current waveform
    % Set up initial guess for shift-scale transform of template
    m = std(cw)/std(template);
    ss0 = [mean(cw)-m*mean(template) m];
    ss(i,:) = fminsearch(@(x) align_shiftscale(x,template,cw),ss0);
    [dist(i),path{i}] = align_shiftscale(ss(i,:),template,cw);  % Get path too
    if (nargin > 2)
      % Calculate x-coordinate of path, in temporal coordinates of the
      % original (unwarped) pitch trajectory
      [utmp,uindx] = unique(path{i}(:,2)); % Find places where pitch indx
                                           % increments
      pitchx{i} = zeros(1,length(uindx));
      for j = 1:length(uindx)-1
        % If several x values are crammed into a single time point,
        % let the returned phase be the average of the values
        pitchx{i}(j) = mean(x(path{i}(uindx(j):uindx(j+1)-1,1)));
      end
      % For the last one, just use the first/only x value
      pitchx{i}(end) = x(path{i}(uindx(end),1));
    end
  end
  