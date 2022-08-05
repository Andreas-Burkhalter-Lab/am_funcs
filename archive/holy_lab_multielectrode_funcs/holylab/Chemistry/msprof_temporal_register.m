function [toffset,crange] = msprof_temporal_register(peakdata,tc,gc,options)
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'N',1000);
  options = default(options,'sigma',options.N/20);
  n_files = length(gc);
  
  %% Put data into matrix format
  % Generate per-sample intensity matrices
  if iscell(peakdata)
    Ic = peakdata;
  elseif isstruct(peakdata)
    n_peaks = length(peakdata);
    Ic = cell(1,n_files);
    for k = 1:n_files
      Ic{k} = zeros(n_peaks,length(peakdata(1).iProf{k}));
    end
    for i = 1:n_peaks
      for k = 1:n_files
        Ic{k}(i,:) = peakdata(i).iProf{k};
      end
    end
  end
  % Determine how to resample in time to get a uniform gradient across
  % files
  [crange,trs] = resample_gradient(gc,options.N);
  % Resample with no shift
  Irs = msptr_resample(Ic,trs,tc,zeros(1,n_files));
  Irs = Irs.^(1/3);  % convert to approx. Gaussian errors

  %% Global registration
  % Register all pairs
  mm = zeros(n_files,n_files);
  for k = 1:n_files
    for kk = k+1:n_files
      [~,p] = register_rigid(Irs(:,:,k),Irs(:,:,kk));
      mm(k,kk) = p(2);
    end
  end
  mm = mm - mm';
  Ioffset = (mean(mm,1) - mean(mm,2)')/2;  % This solution minimizes sum_{i,j} (t_i-t_j-mm(i,j))^2
  dt = mean(diff(trs,1,2),2)';
  toffset = Ioffset.*dt;
  
  %% Local optimization
  % Perform the optimization, i.e., finding the toffset that produces the
  % maximum overlap.
  func = @(toffset) msptr_penalty(Ic,trs,tc,toffset,options.sigma);
  toffset = fminsearch(func,toffset);
  toffset = msptr_shift(toffset,trs,tc);
  
function Irs = msptr_resample(Ic,trs,tc,toffset)
  n_files = length(Ic);
  n_mz = size(Ic{1},1);
  n_t = size(trs,2);
  Irs = zeros(n_mz,n_t-1,n_files);
  % Make sure we don't go over the bounds
  toffset = msptr_shift(toffset,trs,tc);
  % Resample the intensity traces
  for k = 1:n_files
    Irs(:,:,k) = resample_preserve_sum(Ic{k},tc{k},trs(k,:)+toffset(k));
  end
  
function err = msptr_penaltyglobal(t,mm)
  tt = bsxfun(@minus,t(:),t(:)');
  d = tt-mm;
  err = sum(d(:).^2);
  
  
function err = msptr_penalty(Ic,trs,tc,toffset,sigma)
  n_files = length(Ic);
  Irs = msptr_resample(Ic,trs,tc,toffset);
  Irs = Irs.^(1/3);  % convert to approx. Gaussian errors
  % Smooth, so that there is overlap
  Irs = imfilter_gaussian_mex(Irs,[0 sigma 0]);
  % Compute the mismatch. The match is made between individual trace and
  % the sum of all traces, leaving out the "self" trace
  Isum = sum(Irs,3);
  if any(isnan(Isum))
    err = Inf;
    return
  end
  err = 0;
  % The "err" will be -correlation
  for k = 1:n_files
    thisI = Irs(:,:,k);
    Irest = Isum - thisI;
    dp = Irest .* thisI;
    err = err - sum(dp(:))/sqrt(sum(Irest(:).^2)*sum(thisI(:).^2));
  end
  if isnan(err)
    err = Inf;
  end
%   err = err / sum(Irs(:).^2);

function toffset = msptr_shift(toffset,trs,tc)
  % Check to see if any times will be less than 0; if so, shift the rest
  % forward until the smallest is 0
  toffset = toffset - mean(toffset);
  t0rs = trs(:,1) + toffset(:);
  t0 = cellfun(@(p) p(1),tc);
  dtmin = min(t0rs-t0(:))-eps;
  if (dtmin < 0)
    toffset = toffset - dtmin;
  end
  