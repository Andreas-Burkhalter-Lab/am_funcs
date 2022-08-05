function [landmarkClusters, landmarkFinalPos, landmarks, landmarkPos]=clustdata(eventWave, options)
% First step: find a set of landmarks that pretty much cover all the common
% patterns of voltage across electrodes. We use kmeans to do this. After
% kmeans returns, we go through and look for landmarks for which the width
% of the associated points is in the upper half. We then pick another
% landmark from within those "clouds" and re-do kmeans. Thus, the number of
% landmarks increases, but we're biasing their selection so that we get
% more landmarks in areas that are covered more sparsely.

% each n-D point is an event, a scan of all (interested) channels

if(nargin==1)
   options=struct;
end
options=fillOptions(options);

% kmops = struct('showprogress',true);
% [idx,landmarkPos,rmsd] = kmeans_hard(eventWave, options.nInitialLandmarks, kmops); 
% % idx: the landmark ID for each n-D point, that is, column indices to eventWave;
% % landmarkPos: landmark position
% 
% plot(rmsd); drawnow
% 
% n_splits = 2;  % Do the "add new landmarks in sparse areas" twice
% for i = 1:n_splits
%   [landmarks,nlabel] = agglabel(idx);
%   [rmsd,sort_order] = sort(rmsd);
%   landmarkPos = landmarkPos(:,sort_order);
%   landmarks = landmarks(sort_order);
%   nlabel = nlabel(sort_order);
%   ncur = length(rmsd);
%   start = floor(ncur/2);
%   landmarkPos = [landmarkPos nan(size(landmarkPos,1),ncur-start+1)];
%   for j = start:ncur
%     randindx = round(nlabel(j)*rand(1,1)+0.5);
%     landmarkPos(:,ncur+j-start+1) = eventWave(:,landmarks{j}(randindx));
%   end
%   [idx,landmarkPos,rmsd] = kmeans_hard(eventWave,landmarkPos,kmops);
%   plot(rmsd); drawnow
% end

lminfo=choose_landmarks(eventWave, options.nLandmarks); % 
idx=lminfo.landmarkAssignment;
landmarkPos=lminfo.landmarks;

% Now: group together the landmarks that flow uphill to the same place.
% This will organize the landmarks into distinguishable classes.
[ttLandmarkClusters, ttClustinfo] = adaptive_meanshift(eventWave,landmarkPos);
landmarkClusters=agglabel(ttLandmarkClusters);
%      landmarkClusters{i}: all landmarks that form the cluster_i
%      landmarkFinalPos: the settle position for each landmark
landmarkFinalPos=ttClustinfo.yf(:,ttClustinfo.map);

[landmarks,ncl] = agglabel(idx); % Collect all events associated with a
                              % single center/landmark
                              % landmarks{i}: the events(i.e. points) assigned to landmark i
                              
% Compute the average per-electrode s.e.m. associated with each
% landmark. Later we'll use this to pick out electrodes that contain useful
% signal (rather than noise) in identifying the cluster.
clust_rmse_per_electrode = zeros(size(eventWave,1),length(landmarks));
for i = 1:length(landmarks)
  clust_rmse_per_electrode(:,i) = std(eventWave(:,landmarks{i}),1,2)/sqrt(length(landmarks{i}));
end

function options=fillOptions(options)
   if(~isfield(options, 'nLandmarks'))
      options.nLandmarks=5000;
   end
   