function maps = msams_resample(x,varargin)
% Syntax:
%   maps = msams_resample(x)
%   maps = msams_resample(x,y)
%   maps = msams_resample(x,lminfo)
%   maps = msams_resample(...,params)

% Copyright 2008 by Timothy E. Holy

[d,N] = size(x);
if ~isempty(varargin)
  if isnumeric(varargin{1})
    lminfo = choose_landmarks(x,size(varargin{1},2),struct('seed_landmarks',varargin{1}));
  elseif (isstruct(varargin{1}) && ~isfield(varargin{1},'landmarks'))
    lminfo = choose_landmarks(x,N);
  else
    lminfo = varargin{1};
  end
  if (isstruct(varargin{end}) && ~isfield(varargin{end},'landmarks'))
    params = varargin{end};
  else
    params = struct;
  end
end

params = default(params,'resample_size',N,'n_resamples',10);

M = params.resample_size;
n = params.n_resamples;
nLandmarks = size(lminfo.landmarks,2);

maps = nan(n,nLandmarks);
for i = 1:n
  if (M < N/2)
    % Sample without replacement
    ri = randperm(N);
    ri = ri(1:M);
  else
    % Sample with replacement
    ri = round(N*rand(1,M)+0.5);
  end
  x_rs = x(:,ri);
  lminfo_rs = reindex_landmarks(lminfo,ri);
  [clust,clustinfo_rs] = msams(x_rs,lminfo_rs);
  maps(i,:) = clustinfo_rs.map;
end

%     if options.use_mex && have_lminfo
%       clustinfo2 = msams_converge_mex(x(:,ri),lminfo,yf(:,maps_to_self),options);
%     else
%       [yftmp,maptmp,ntmp,clustinfo2] = msams_converge(x,yf(:,maps_to_self),options);
%     end


%   %
%   % Make sure the final assignments are stable by doing a bootstrap
%   % resampling of the underlying distribution
%   %
%   n_self_Old = Inf;
%   maps_to_self = find(map == 1:length(map));
%   n_self = length(maps_to_self);
%   while n_self ~= n_self_Old
%     n_self_Old = n_self;
%     ri = round(N*rand(1,N)+0.5);
%     ri(ri > N) = N;  % for safety, but probably not necessary
%     if options.use_mex && have_lminfo
%       clustinfo2 = msams_converge_mex(x(:,ri),lminfo,yf(:,maps_to_self),options);
%     else
%       [yftmp,maptmp,ntmp,clustinfo2] = msams_converge(x,yf(:,maps_to_self),options);
%     end
%     mapOld = map;
%     map(maps_to_self) = clustinfo2.closestLandmark;
%     map = leapfrog(map);
%     maps_to_self = find(map == 1:length(map));
%     n_self = length(maps_to_self);
%   end
%   clustinfo.map = map;
%   % Re-label the clusters 1,2,3,...
%   [umap,tmp,clust] = unique(map);
%   return