function [componentsm,lambda] = components_from_waveform(v,options)
  [nchan,N] = size(v);
  if (nargin < 2)
    options = struct;
  end
  required_fieldnames = {'sniprange','polarity','n_components'};
  missing_fieldnames = setdiff(required_fieldnames,fieldnames(options));
  if ~isempty(missing_fieldnames)
    fprintf('Missing fields:\n');
    for i = 1:length(missing_fieldnames)
      fprintf('%s\n',missing_fieldnames{i});
    end
    error('Must have all required fieldnames');
  end
  options = default(options,'n_components',round(1.5*nchan),'max_snips',Inf,'isMergeAdjacent',true,'adjacentNScans',5,'maxIters',2,'thresh_factor',6,'skip_baselineshift',false);
  
  % Compute some statistics on the input waveform
  if ~options.skip_baselineshift
    vmed = median(v,2);
    v = v - repmat(vmed,1,N);
  end
  vnoise = mean(abs(v),2); % robust estimation of noise
  
  if ~isfield(options,'thresh')
    % Compute thresholds
    options.thresh = options.thresh_factor*vnoise;
  end
  
  % Find spike times
  fprintf('Finding spike times...');
  t = findpeaks_multichan(v,options.thresh,options);
  % Toss any times that are too close to edges
  t = snipedges(t,N,options);
  % Choose a subset of times, if there are too many
  if (length(t) > options.max_snips)
    skip = ceil(options.max_snips/length(t));
    t = t(1:skip:end);
  end
  n_snips = length(t);
  sniplen = diff(options.sniprange)+1;
  rng = options.sniprange(1):options.sniprange(2);
  rng0 = 0:sniplen-1;
  fprintf('done\n');
  
  % Scale each channel by its "noise"
  % (This is a "cheap" method for normalizing by the noise that does not
  % run into problems of limited statistics that might come from a real
  % estimate of the noise covariance matrix.)
  %vrms = std(v,0,2);
  vnorm = vnoise;
  for i = 1:nchan
    v(i,:) = v(i,:) / vnorm(i);
  end

  % Check snippets to make sure none have NaNs (which would be used to
  % encode undefined regions/boundaries)
  badFlag = false(1,n_snips);
  for i = 1:n_snips
    snip = v(:,t(i) + rng);
    if any(isnan(snip(:)))
      badFlag(i) = true;
    end
  end
  
  t(badFlag) = [];
  n_snips = length(t);
  tshift = t + options.sniprange(1);
  fprintf('%d spike snippets were cut.\n',n_snips);
  
  % Compute the covariance matrix of these events
  fprintf('Computing covariance...');
  C = covar_from_waveform(v,t,options.sniprange);
  fprintf('done\n');
  
  % In preparation for the overlap decomposition, split the times into
  % non-overlapping blocks (will reduce the memory requirements).
  fprintf('Separating into non-overlapping blocks...');
  dt = diff(t);
  splitIndex = [0 find(dt > sniplen) length(t)];
  % Don't do _too_ many separate blocks, it will be inefficient
  ilast = 1;
  inext = 3;
  killFlag = false(1,length(splitIndex));
  while (inext <= length(splitIndex))
     while (inext <= length(splitIndex) && splitIndex(inext) - splitIndex(ilast) < 100)
       killFlag(inext-1) = true;
       inext = inext+1;
     end
     ilast = inext;
     inext = inext+1;
  end
  splitIndex(killFlag) = [];
  fprintf('done (%d blocks used)\n',length(splitIndex)-1);
  
  % Calculate the components. We do this iteratively, allowing the
  % improvements in the ability to subtract off overlapping spikes to
  % inform the construction of the components.
  for iter = 1:options.maxIters
    % Do the eigenvalue decomposition
    fprintf('Doing the eigenvalue decomposition...');
    [components,D] = eig(C);
    fprintf('done\n');
    lambda = diag(D);
    lambda = lambda(end:-1:1);
    if (iter < options.maxIters)
      n_components = min([n_snips options.n_components nchan*sniplen]);
    else
      % Keep them all on the final iteration (so user can change mind about
      % how many to keep!)
      n_components = size(components,2);
    end
    components = components(:,end:-1:end-n_components+1);
    componentsm = vec2snip(components,nchan);
    if (iter < options.maxIters)
      % Calculate the amplitude of each component for each spike. We don't
      % use the SVD result here, because we want to do it in a way that
      % properly handles overlaps.
      fprintf('Computing component overlaps...');
      c_overlap = component_overlaps(componentsm);
      fprintf('done\n');
      fprintf('Computing the right hand side of the linear equation...');
      b = wdecomp_amplitude_rhs(v,componentsm,tshift);
      fprintf('done\n');
      a = zeros(size(b));
      % Split into non-overlapping chunks, to avoid using too much memory
      fprintf('Computing the hessian and solving the linear equation')
      for i = 1:length(splitIndex)-1
        istart = splitIndex(i)+1;
        iend = splitIndex(i+1);
        H = wdecomp_amplitude_hessian(c_overlap,tshift(istart:iend));
        btmp = b(:,istart:iend);
        atmp = H \ btmp(:);
        a(:,istart:iend) = reshape(atmp,[n_components iend-istart+1]);
        fprintf('.');
      end
      fprintf('done\n');
      % Calculate the residual.  This will allow us to (approximately)
      % add back in the "noise" (to try to give "real" snippets) while
      % still largely eliminating overlapping spikes.  This should give
      % waveforms for the next SVD that will be largely uncontaminated by
      % overlaps.
      fprintf('Calculating the residual...');
      res = v;
      for eventIndex = 1:n_snips
        % Calculate the shape of the waveform snippet from the components
        % and their amplitudes---this gives an estimate of the "true" shape
        % without overlapping events.
        snip = components * a(:,eventIndex);
        snip = reshape(snip,[nchan sniplen]);
        % Subtract snippet from waveform to get residual
        res(:,tshift(eventIndex) + rng0) = res(:,tshift(eventIndex) + rng0) - snip;
      end
      fprintf('done\n');
      % Compute the snippet covariance matrix, correcting each
      % component-based snippet by the residual.
      % (in cases where a snippet is isolated from others, this just
      % "undoes" the residual calculation, but for cases with overlap
      % going through the residual was an important step)
      fprintf('Recalculating the covariance of overlap-corrected spike waveforms...');
      C = covar_from_waveform(res,t,options.sniprange,double(components),double(a));
      fprintf('done\n');
    end
  end
  %components = components .* vnorm(:,ones(1,n_components));
  %components = components ./ repmat(sqrt(sum(components.^2,1)),[nchan*sniplen 1]);
  % Scale each channel back to original size
  componentsm = componentsm .* repmat(vnorm,[1 sniplen n_components]);
  % Normalize
  cnorm = sqrt(sum(sum(componentsm.^2,1),2));
  componentsm = componentsm ./ repmat(cnorm,[nchan sniplen 1]);
  
    
  
function t = snipedges(t,N,options)
  t = t(t > max(0,-options.sniprange(1)));
  t = t(t < N - max(0,options.sniprange(2)));

% function vec = snip2vec(snip)
%   sz = size(snip);
%   vec = reshape(snip,[sz(1)*sz(2) sz(3)]);
  
function snip = vec2snip(vec,nchannels)
  [nprod,nsnips] = size(vec);
  snip = reshape(vec,[nchannels nprod/nchannels nsnips]);
  