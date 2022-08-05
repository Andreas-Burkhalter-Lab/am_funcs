function funcs = fit_utilities
% FIT_UTILITIES: a set of convenience functions for spike sorting
% Note: these are poorly documented; look at the source code. Their main
% purpose is to be used internally in spike sorting code and for
% debugging. Also, a warning: the API for these functions could
% change.
%
% Usage:
%   funcs = fit_utilities
% will return a structure with a set of function handles to the
% convenience functions. (This prevents the need for a whole lot of
% tiny functions cluttering up the directory.)
  
% Copyright 2008 by Timothy E. Holy
  
  funcs.flag_times_in_range = @fu_flag_times_in_range;

  funcs.snip2vec_by_t = @fu_snip2vec_by_t;
  funcs.snip2vec_by_c = @fu_snip2vec_by_c;
  funcs.vec2snip_by_t = @fu_vec2snip_by_t;
  funcs.vec2snip_by_c = @fu_vec2snip_by_c;
  
  funcs.plot_snips_on_array = @fu_plot_snips_on_array;
  
  funcs.collect_snippets = @fu_collect_snippets;

  funcs.fit_waveform = @fu_fit_waveform;
  funcs.solve_for_amplitudes = @fu_solve_for_amplitudes;

  funcs.calculate_residual = @fu_calculate_residual;
  funcs.reconstruct_snippets_using_residual = @fu_reconstruct_snippets_using_residual;
  funcs.minmax_using_residual = @fu_minmax_using_residual;
  funcs.residual_errors = @fu_residual_errors;
  
  funcs.run_clustering = @fu_run_clustering;
  
  funcs.svd_components = @fu_svd_components;
  funcs.split_times = @fu_split_times;
end


% General notes: for much of what follows, a "snip" is a
% nchans-by-sniplen-by-nsnips array.

function keepFlag = fu_flag_times_in_range(t,range,sniprange)
  if (nargin > 2)
    keepFlag = t + sniprange(1) >= range(1) & ...
	t + sniprange(2) <= range(2);
  else
    keepFlag = t >= range(1) & t <= range(2);
  end
end

% The "by_c" variants of the following functions are much faster (they
% don't have an extra "permute" step), but are less convenient for
% plotting. So, the "by_t" variants are especially useful while debugging.

function vec = fu_snip2vec_by_c(snip)
  sz = [size(snip) 1];
  vec = reshape(snip,[sz(1)*sz(2) sz(3)]);
end

function vec = fu_snip2vec_by_t(snip)
  sz = [size(snip) 1];
  snip = permute(snip,[2 1 3]);
  vec = reshape(snip,[sz(1)*sz(2) sz(3)]);
end

function snip = fu_vec2snip_by_c(vec,nchan)
  sz = size(vec);
  if (sz(1)/nchan ~= round(sz(1)/nchan))
    error('Size mismatch');
  end
  snip = reshape(vec,[nchan sz(1)/nchan sz(2)]);
end

function snip = fu_vec2snip_by_t(vec,nchan)
  sz = size(vec);
  if (sz(1)/nchan ~= round(sz(1)/nchan))
    error('Size mismatch');
  end
  snip = reshape(vec,[sz(1)/nchan nchan sz(2)]);
  snip = permute(snip,[2 1 3]);
end

% An especially convenient syntax for plotting spikes in the array geometry

function hax_out = fu_plot_snips_on_array(snip,arrayname,hax,options)
  if (nargin < 3)
    figure;
    hax = gca;
  end
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'showscale',true);
  snip = permute(snip,[2 1 3]);
  ch = get_hda_holylab(arrayname);
  snip = shape_spikes_to_array(snip,ch,ch);
  [hline,hax_outo] = plot_spikes_on_array(snip,hax,options);
  if (nargout > 0)
    hax_out = hax_outo;
  end
end


function snip = fu_collect_snippets(v,t,sniprange)
  nchan = size(v,1);
  n_snips = length(t);
  sniplen = diff(sniprange)+1;
  rng = sniprange(1):sniprange(2);

  snip = zeros(nchan,sniplen,n_snips);
  for snipIndex = 1:n_snips
    thisT = t(snipIndex);
    snip(:,:,snipIndex) = v(:,thisT + rng);
  end
end  

% Don't use the following repetitively; it's inefficient because of the
% call to component_overlaps. This is useful mostly for one-shot deals
function a = fu_fit_waveform(v,tshift,componentsAsSnips,progress)
  % tshift is tpeak + sniprange(1)
  if (nargin < 4)
    progress = false;
  end
  if progress
    fprintf('Computing component overlaps...');
  end
  overlap = component_overlaps(componentsAsSnips);
  if progress
    fprintf('done\n');
    fprintf('Computing the right hand side of the linear equation...');
  end
  b = wdecomp_amplitude_rhs(v,componentsAsSnips,tshift);
  if progress
    fprintf('done\n');
  end
  a = fu_solve_for_amplitudes(tshift,b,overlap,100,progress);
end

% A version that allows you to reduce memory requirements, and also _can_
% be called repetitively because you compute the overlap matrix externally
function a = fu_solve_for_amplitudes(tshift,b,overlap,nmax,progress)
  % tshift is tpeak + sniprange(1)
  if (nargin < 5)
    progress = false;
  end
  if (nargin < 4)
    nmax = Inf;
  end
  [n_components n_spikes] = size(b);
  if (length(tshift) ~= n_spikes)
    error('Size mismatch between b and t');
  end
  if (size(overlap,1) ~= n_components)
    error('Size mismatch between overlap and b');
  end
  if (length(tshift) > nmax)
    % Break times into non-overlapping blocks to reduce memory requirements
    sniplen = (size(overlap,3)+1)/2;
    [splitRange,saveRange] = fu_split_times(tshift,sniplen,nmax);
    if progress
      fprintf('Computing the hessian and solving the linear equation. %% done: ')
    end
    a = zeros(size(b));
    H1 = full(wdecomp_amplitude_hessian(overlap,0));
    pctdone_last = 0;
    for i = 1:size(splitRange,2)
      istart = splitRange(1,i);
      iend = splitRange(2,i);
      if (istart == iend)
        a(:,istart) = H1\b(:,istart);
      else
        H = wdecomp_amplitude_hessian(overlap,tshift(istart:iend));
        btmp = b(:,istart:iend);
        atmp = H \ btmp(:);
        atmp = reshape(atmp,[n_components iend-istart+1]);
        rngsave = saveRange(1,i):saveRange(2,i);
        keepIndex = findainb(rngsave,istart:iend);
        atmp = atmp(:,keepIndex);
        a(:,rngsave) = atmp;
      end
      if progress
        pctdone = floor(istart/splitRange(end)*100);
        if pctdone > pctdone_last
          fprintf('...%d',pctdone);
          pctdone_last = pctdone;
        end
      end
    end
    if progress
      fprintf('done\n');
    end
  else
    % This is the path for "short" stretches that can be done in a single go
    if progress
      fprintf('Calculating the Hessian...');
    end
    H = wdecomp_amplitude_hessian(overlap,tshift);
    if progress
      fprintf('done\nSolving the linear equation...');
    end
    a = H \ b(:);
    a = reshape(a,[n_components n_spikes]);
    if progress
      fprintf('done\n');
    end
  end
end

function [nchan,sniplen,n_components,n_snips] = fu_validate_component_inputs(v,tshift,componentsAsSnips,a)
  [nchan,sniplen,n_components] = size(componentsAsSnips);
  n_snips = length(tshift);
  if (size(v,1) ~= nchan)
    error('Size mismatch between residual and components');
  end
  if (nargin > 3)
    if ~isequal(size(a),[n_components,n_snips])
      error('Size mismatch between a and other inputs');
    end
  end
end

function res = fu_calculate_residual(v,tshift,componentsAsSnips,a)
  % tshift is tpeak + sniprange(1)
  [nchan,sniplen,n_components,n_snips] = fu_validate_component_inputs(v,tshift,componentsAsSnips,a);
  
  res = v;
  componentsv = fu_snip2vec_by_c(componentsAsSnips);
  rng0 = 0:sniplen-1;
  for eventIndex = 1:n_snips
    % Calculate the shape of the waveform snippet from the components
    % and their amplitudes---this gives an estimate of the "true" shape
    % without overlapping events.
    snip = componentsv * a(:,eventIndex);
    snip = reshape(snip,[nchan sniplen]);
    % Subtract snippet from waveform to get residual
    res(:,tshift(eventIndex) + rng0) = res(:,tshift(eventIndex) + rng0) - snip;
  end
end

function snip = fu_reconstruct_snippets_using_residual(res,tshift,componentsAsSnips,a)
  % tshift is tpeak + sniprange(1)
  [nchan,sniplen,n_components,n_snips] = fu_validate_component_inputs(res,tshift,componentsAsSnips,a);

  componentsv = fu_snip2vec_by_c(componentsAsSnips);
  rng0 = 0:sniplen-1;
  snip = zeros([nchan sniplen n_snips],'single');
  for eventIndex = 1:n_snips
    tsnip = componentsv * a(:,eventIndex);
    tsnip = reshape(tsnip,[nchan sniplen]);
    snip(:,:,eventIndex) = res(:,tshift(eventIndex) + rng0) + tsnip;
  end
end

% This also gives the same output as residual_errors, if desired
function [pp,err] = fu_minmax_using_residual(res,tshift,componentsAsSnips,a)
  [nchan,sniplen,n_components,n_snips] = fu_validate_component_inputs(res,tshift,componentsAsSnips,a);
  
  componentsv = fu_snip2vec_by_c(componentsAsSnips);
  rng0 = 0:sniplen-1;
  pp = zeros([nchan 2 n_snips],'single');
  err = zeros(1,n_snips);
  for eventIndex = 1:n_snips
    thisres = res(:,tshift(eventIndex) + rng0);
    tsnip = componentsv * a(:,eventIndex);
    tsnip = reshape(tsnip,[nchan sniplen]);
    tsnip = thisres + tsnip;
    pp(:,1,eventIndex) = min(tsnip,[],2);
    pp(:,2,eventIndex) = max(tsnip,[],2);
    if (nargout > 1)
      err(eventIndex) = sum(thisres(:).^2);
    end
  end
  if (nargout > 1)
    err = sqrt(err);
  end
end

function err = fu_residual_errors(res,t,sniprange)
  rng = sniprange(1):sniprange(2);
  n_snips = length(t);
  err = zeros(1,n_snips);
  for eventIndex = 1:n_snips
    tmp = sum(res(:,t(eventIndex)+rng).^2);
    err(eventIndex) = sum(tmp);
  end
  err = sqrt(err);
end

function clust = fu_run_clustering(x,maxlm,msamsops)
  if (nargin < 2)
    maxlm = 5000;
  end
%   % First find the "true" landmarks that will speed clustering
%   nlm = min(maxlm,ceil(sqrt(size(x,2))));
%   lminfo = choose_landmarks(x,nlm);
%   % Now choose the remaining probe points using the same algorithm
%   lminfo2 = choose_landmarks(x,maxlm,lminfo);
  lminfo = choose_landmarks(x,maxlm);
  if (nargin < 3)
    msamsops = struct;
  end
  msamsops = default(msamsops,'reflow',true,'consolidate',false,'flow_only_landmarks',true);
  clust = msams(x,lminfo,msamsops);
  % Apply to all data points, not just the probe points
%   clust = clust(lminfo2.landmarkAssignment);
end

function [csvd,s] = fu_svd_components(componentsAsSnips)
  cv = fu_snip2vec_by_c(componentsAsSnips);
  [U,S] = svd(cv,'econ');
  s = diag(S);
  csvd = fu_vec2snip_by_c(U,size(componentsAsSnips,1));
end

function [splitRange,saveRange] = fu_split_times(tshift,sniplen,nmax)
  % Splits times up in to chunks (perhaps with overlap) for which the
  % linear equation can be productively solved without overloading memory
  n_spikes = length(tshift);
  dt = diff(tshift);
  splitIndex = [0 find(dt > sniplen) length(tshift)]; % events separated by more than sniplen are independent
  
  % Consolidate tiny blocks, to avoid being too inefficient
  ilast = 1;
  inext = 3;
  killFlag = false(1,length(splitIndex));
  while (inext <= length(splitIndex))
    while (inext <= length(splitIndex) && splitIndex(inext) - splitIndex(ilast) < nmax)
      killFlag(inext-1) = true;
      inext = inext+1;
    end
    ilast = inext-1;
    inext = inext+1;
  end
  %splitIndex(killFlag) = [];
  splitRange = [splitIndex(1:end-1)+1; splitIndex(2:end)];
  saveRange = splitRange;
  
  % In cases where a block is of size larger than nmax, we have to
  % break it. We break it with some overlap, so that we don't get
  % significant inaccuracies at the edges.
  tooBigIndex = find(diff(splitRange,1,1) > nmax);
  for i = tooBigIndex(end:-1:1)
    newIndex = splitRange(1,i)-1:nmax:splitRange(2,i);
    if (newIndex(end) ~= splitRange(2,i))
      newIndex = [newIndex splitRange(2,i)];
    end
    newRangeSave = [newIndex(1:end-1)+1; newIndex(2:end)];
    newRangeCalc = newRangeSave;
    for j = 1:size(newRangeSave,2)-1
      % Shift the right edge of the calculation region farther to the
      % right. Do this so that at least 3*sniplen samples are used as the
      % "edge"
      isok = false;
      while ~isok
        isok = tshift(newRangeCalc(2,j)) - tshift(newRangeSave(2,j)) > 3*sniplen;
        if ~isok
          if newRangeCalc(2,j) < n_spikes
            newRangeCalc(2,j) = newRangeCalc(2,j)+1;
          else
            isok = true;
          end
        end
      end
      % Shift the left edge of the calculation region farther to the left
      isok = false;
      while ~isok
        isok = tshift(newRangeSave(1,j+1)) - tshift(newRangeCalc(1,j+1)) > 3*sniplen;
        if ~isok
          if newRangeCalc(1,j+1) > 1
            newRangeCalc(1,j+1) = newRangeCalc(1,j+1)-1;
          else
            isok = true;
          end
        end
      end
    end
    splitRange = [splitRange(:,1:i-1) newRangeCalc splitRange(:,i+1:end)];
    saveRange = [saveRange(:,1:i-1) newRangeSave saveRange(:,i+1:end)];
  end
end

