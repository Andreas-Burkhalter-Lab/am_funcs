function [Ms,minmax,mapf,selfIndex,domains] = mspeaks_by_imflow(varargin)
  calculateM = true;
  calculate_map = true;
  if isstruct(varargin{1})
    scan = varargin{1};
    min_mz = min([scan.min_mz]);
    max_mz = max([scan.max_mz]);
    n_bins_mz = round((max_mz - min_mz)*4);
    if (length(varargin) < 2)
      options = struct;
    else
      options = varargin{2};
    end
  else
    Ms = varargin{1};
    n_bins_mz = Ms;
    calculateM = false;
    minmax = [1 size(Ms,1)];
    carg = 2;
    options = struct;
    while (carg <= length(varargin))
      if isnumeric(varargin{carg})
        mapf = varargin{carg};
        calculate_map = false;
      elseif isstruct(varargin{carg})
        options = varargin{carg};
      end
      carg = carg+1;
    end
  end
  options = default(options,'n_bins_mz',n_bins_mz,...
    'smoothing',[0.5 20],...
    'thresh',1e4);
  
  if calculateM
    fprintf('Extracting & binning the data...');
    [M,minmax] = msdecimate(scan,struct('n_bins_mz',options.n_bins_mz));
    fprintf('done.\nSmoothing...');
    Ms = imfilter_gaussian(M,options.smoothing);
    fprintf('done.\n');
    clear M  % We have to clear as we go for memory reasons
  end
  
  if calculate_map
    % Find the flow map
    fprintf('Finding the map...');
    map = imflow(Ms);
    fprintf('done.\n')
    % Flow the map until the peaks are identified
    fprintf('Flowing the map...');
    mapf = map;
    counter = 0;
    mapfOld = 0*mapf;
    while ~isequal(mapf,mapfOld)
      mapfOld = mapf;
      mapf = mapf(mapf);
      counter=counter+1;
    end
    fprintf('done (required %d iterations).\n',counter);
    mapf = killedges(mapf,0);
    clear map mapfOld
  end

  % Find the points at the peak
  self = (mapf == reshape(1:numel(mapf),size(mapf)));
  
  % Find the big peaks
  selfIndex = find(self);
  v = Ms(selfIndex);
  selfIndex = selfIndex(v > options.thresh);
  
  % Find the domains
  fprintf('Finding the domains...');
  [umap,tmp,domainIndex] = unique(mapf(:));
  domains = agglabel(domainIndex);
  [ii,dIndx] = intersect(umap,selfIndex);  % keep domains that are in selfIndex
  domains = domains(dIndx);
  fprintf('done.\n');
end