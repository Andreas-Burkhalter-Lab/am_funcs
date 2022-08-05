function [gain,bias,noise] = imnoiseanalysis(ip,options)
% IMNOISEANALYSIS: analyze properties of cameras and PMTs
% A series of images (with repeated measurement of one or more
% nominally-identical scenes) are input at different gains and
% integration times.  These images then allow one to calculate the camera
% gain (conversion factor between electrons and digital numbers), bias,
% and noise.
% 
% This function employs the following model of pixel intensity p:
%     p = g(G) * nsignal + bias + noise
% where the bias term is constant (for each pixel, gain combination) and
% the noise is gaussian with zero mean (for each pixel, gain
% combination). g(G) is the _actual_ gain of the system (electrons to
% digital numbers), as a function of the user-supplied gain setting G.
% nsignal is rsignal*t, where t is the integration time; rsignal includes
% photon-emission by the sample (if any) and dark current. nsignal is
% assumed to be Poisson, i.e., the variance is equal to the mean.
%
% Guide to image acquisition:
%   1. If all you want is the gain, then it's not necessary to vary any
%      of the acquisition parameters, as long as you're confident that the
%      bias is uniform across the chip: just take some (minimum 2)
%      nominally-identical images of a scene containing regions of
%      varying intensity.
%   2. Bias and signal can be separated only by taking images with varying
%      exposure times, so do that if you want to extract these
%      parameters.  Minimally, you need 2 different exposure times.
%      In general, you want to avoid having a lot of saturated pixels.
%   3. If you supply images with multiple gain settings, you can extract
%      these parameters as a function of gain setting.
%
% Syntax:
%   [gain,bias,noise] = imnoiseanalysis(ip,options)
%   [....] = imnoiseanalysis(ip,options)  
% where
%   ip is an input IMPHYS structure.  If all the settings are the same
%     (so that a set of nominally-identical frames are being supplied),
%     then no additional information needs to be supplied.  Only the gain
%     will be computed.
%     If you want to supply frames taken with different settings, you
%     specify the settings by adding the following fields to the IMPHYS
%     structure:
%       gainsetting: the user-supplied gain setting
%       inttime: the integration time
%     If you omit gainsetting, it assumes all are taken with the same
%     gainsetting.
%   options is an optional structure with the following fields:
%     plotbias (default 500): if > 0, randomly selects plotbias pixels and
%       plots their mean pixel value vs. integration time;
%     plotgain (default true): if true, shows the scatter plot of
%       var(pixel) vs. mean(pixel);
%     numbins (default 100): sets the number of bins used in the
%       variance-to-mean plot;
%     discardmax (default true): if true, throws out any pixel which
%       reaches the maximum value returned by the camera (saturated
%       pixels), as reported by the imrange parameter;
%     discardmin (default false): if true, throws out any pixel which
%       reaches the minimum value returned by the camera
%     meanfrac (default [0.1 0.9]): allows you to include only a subset
%       of the pixels in doing the regression---the 2-vector yields the
%       percentile values, so [0.5 1] includes the top-half of the pixel
%       values in the regression.  Set to [0 1] to use all pixels.
%
% The outputs are as follows:
%   gain is a scalar yielding the gain, if only one gain setting is
%     supplied. If there are multiple gain settings, this is a
%     nsettings-by-2 matrix, where the first column is the gain setting,
%     and the second column is the actual gain.
%     different gain setting.
%   bias is a matrix of the dimensions of a single image, or a stack of
%     frames if different gain settings are provided (one for each
%     setting).
%   noise is the mean variance of the noise (one for each gain setting).
%
% See also: ANALYZECAMERA
  
% Copyright 2004 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  options = imnaoptions(options);

  % Get image dimensions
  nx = diff(ip(1).xrange)+1;
  ny = diff(ip(1).yrange)+1;
  npix = nx*ny;
  % Collect frames with particular gain setting & integration time
  if (isfield(ip,'gainsetting') & isfield(ip,'inttime'))
    gaintime = [[ip.gainsetting]; [ip.inttime]]';
  elseif isfield(ip,'gainsetting')
    gaintime = [[ip.gainsetting]', ones(length(ip),1)];
  elseif isfield(ip,'inttime')
    gaintime = [ones(length(ip),1), [ip.inttime]'];
  else
    gaintime = ones(length(ip),2);
  end
  [ugaintime,uindx,gaintimeindx] = unique(gaintime,'rows');
  
  % Determine whether multiple integration times were supplied
  % (and send error if not and user is asking for too much)
  [ut,uindx,ugttindx] = unique(ugaintime(:,2));
  if (length(ut) == 1 & nargout > 1)
    error(['You must supply data from multiple integration times if you ' ...
           'want bias and noise information']);
  end

  % Compute means and variances of frames having identical settings
  nu = size(ugaintime,1);   % Number of unique settings
  fmean = zeros(ny,nx,nu);
  fvar = zeros(ny,nx,nu);
  flagbad = zeros(ny,nx);
  for i = 1:nu
    fprintf('Settings %d of %d\n',i,nu);
    indx = find(gaintimeindx == i);  % These are the frames with current
                                     % settings
    nreps = length(indx);
    if (nreps > 0)
      tmp = imphysfetch(ip(indx(1)));
      frametmp = zeros(ny,nx,nreps,class(tmp));
      frametmp(:,:,1) = tmp;
      for j = 2:nreps
        frametmp(:,:,j) = imphysfetch(ip(indx(j)));
      end
      [fvar(:,:,i),fmean(:,:,i)] = varmean(frametmp,0,3);
    end
    for j = 1:nreps
      if options.discardmax
        indxbad = find(frametmp(:,:,j) >= ip(indx(1)).imrange(2));
        flagbad(indxbad) = 1;
      end
      if options.discardmin
        indxbad = find(frametmp(:,:,j) <= ip(indx(1)).imrange(1));
        flagbad(indxbad) = 1;
      end
    end
  end
  % Reshape mean & var to "eliminate" spatial info
  fmean = reshape(fmean,[npix nu])';
  fvar = reshape(fvar,[npix nu])';
  flagbad = reshape(flagbad,[npix 1])';
  indxgood = find(flagbad == 0);
  if ~isempty(find(flagbad))
     warning([num2str(length(find(flagbad))) ' saturated pixels found']);
  end
  
  % There are two cases to be treated: one in which the bias can be
  % computed (because there are at least 2 integration times), and one in
  % which it cannot.  Handle the "cannot" case first; wing it on the
  % gain computation, assuming the bias is uniform.
  [ug,uindx,ugttindx] = unique(ugaintime(:,1));  % Combine frames with
                                                 % same gain setting
  gain = [ug,ug];
  ngains = length(ug);
  if (length(ut) == 1)
    options.biassubtract = 0;
    for i = 1:ngains
      tindx = find(ugttindx == i);
      mtmp = fmean(tindx,indxgood);
      vtmp = fvar(tindx,indxgood);
      % gain = slope of variance vs. mean line
      gain(i,2) = imna_gain_noise(mtmp,vtmp,options);
    end
    if (ngains == 1)
      gain = gain(2);
    end
    return
  end
  
  
  % OK, let's do the full computation
  options.biassubtract = 1;
  biasgood = zeros(ngains,length(indxgood));
  noise = zeros(1,ngains);
  for i = 1:ngains   % Compute for each gain separately
    tindx = find(ugttindx == i);
    mtmp = fmean(tindx,indxgood);
    % bias = mean pixel value extrapolated to 0 integration time
    if options.plotbias
      figure
      rp = randperm(size(mtmp,2));
      %loglog(ugaintime(tindx,2),mtmp(:,rp(1:options.plotbias)));
      plot(ugaintime(tindx,2),mtmp(:,rp(1:options.plotbias)));
      xlabel('Integration time')
      ylabel('Mean pixel value')
      title(['Gain setting ' num2str(gain(i,1))]);
    end
    [slope,biasgood(i,:)] = linregress(ugaintime(tindx,2),mtmp);
    mtmp = mtmp - repmat(biasgood(i,:),length(tindx),1); % bias-subtracted vals
    vtmp = fvar(tindx,indxgood);
    [gain(i,2),noise(i)] = imna_gain_noise(mtmp,vtmp,options,gain(i,1));
    % noise = var - gain*(bias-subtracted mean value)
    %noise(i,:) = mean(vtmp - gain(i,2)*mtmp);
  end
  
  % Reshape the outputs
  if (ngains == 1)
    gain = gain(2);
    %bias = reshape(bias,[ny nx]);
    bias = nan(ny,nx);
    bias(indxgood) = biasgood;
  else
    %bias = reshape(bias',[ny nx ngains]);
    bias = nan(npix,ngains);
    for i = 1:ngains
       bias(indxgood,i) = biasgood(i,:)';
    end
    bias = reshape(bias,[ny nx ngains]);
  end
    
  
function options = imnaoptions(options)
  if ~isfield(options,'plotbias')
    options.plotbias = 500;
  end
  if ~isfield(options,'plotgain')
    options.plotgain = 1;
  end
  if ~isfield(options,'meanfrac')
    options.meanfrac = [0.1 0.9];
  end
  if ~isfield(options,'discardmax')
    options.discardmax = 1;
  end
  if ~isfield(options,'discardmin')
    options.discardmin = 0;
  end
  if ~isfield(options,'numbins')
    options.numbins = 100;
  end
  
  
function [g,n] = imna_gain_noise(mtmp,vtmp,options,gainset)
  % gain = slope of variance vs. mean line
  % noise = y-intercept of variance vs. mean line
  [smtmp,sortindx] = sort(mtmp(:));
  regressrng = round(options.meanfrac*(length(sortindx)-1))+1;
  [g,n] = linregress(smtmp(regressrng(1):regressrng(2)), ...
                     vtmp(sortindx(regressrng(1):regressrng(2))));
  if options.plotgain
    % To make the figure manageable, do some averaging to reduce the
    % number of points
    splitpoints = round(linspace(1,length(smtmp),options.numbins));
    svtmp = vtmp(sortindx);
    for k = 1:length(splitpoints)-1
      px(k) = mean(smtmp(splitpoints(k):splitpoints(k+1)));
      py(k) = mean(svtmp(splitpoints(k):splitpoints(k+1)));
      pye(k) = std(svtmp(splitpoints(k):splitpoints(k+1)))/ ...
        sqrt(diff(splitpoints([k k+1])));
    end
    figure
    %h = errorbar(px,py,pye);
    %set(h,'LineStyle','none');
    h = line([px; px],[py-pye;py+pye],'Color','b','LineWidth',2);
    if options.biassubtract
      xlabel('Bias-subtracted mean pixel value')
    else
      xlabel('Mean pixel value')
    end
    ylabel('Var pixel value')
    pym = px * g + n;
    line(px,pym,'Color','r');
    if ~options.biassubtract
      n = NaN;
    end
    if (nargin > 3)
      title(sprintf('Gain setting %g: gain=%g, var(noise)=%g => %g e- r.m.s',...
        gainset,g,n,sqrt(n/g^2)));
    else
      title(sprintf('Gain=%g, var(noise)=%g',g,n));
    end
  end
    
