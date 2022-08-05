function [gain,bias,noise] = imnoiseanalysis(smm,varargin)
% IMNOISEANALYSIS: analyze properties of cameras and PMTs
% A series of images (with repeated measurement of one or more
% nominally-identical scenes) are input at different gains and
% integration times.  These images then allow one to calculate the camera
% gain (conversion factor between electrons and digital numbers), bias,
% and noise.
% 
% This function employs the following model of pixel intensity p:
%     p = g * nsignal + bias + noise
% where the bias term is constant (for each pixel, gain combination) and
% the noise is gaussian with zero mean (for each pixel, gain
% combination). "g" is the _actual_ gain of the system (converts quantal
% units to gray levels, so it has units of grays/quantum), which of course
% depends upon the user settings. "nsignal" is rsignal*t, where t is the
% integration time; rsignal includes photon-emission by the sample (if any)
% and dark current. nsignal is assumed to be Poisson, i.e., the variance is
% equal to the mean. "noise" is assumed to be Gaussian-distributed, and has
% units of grays^2.  "bias" may either be constant across the image, or be defined
% on a pixel-by-pixel basis, and has units of grays.
%
% Guide to image acquisition:
%   1. If all you want is the gain, then it's not necessary to vary any
%      of the acquisition parameters, as long as you're confident that the
%      bias is uniform across the chip: just take some (minimum 2,
%      preferably more) nominally-identical images of a scene containing
%      regions of varying intensity.
%   2. Bias and signal can be separated only by taking images with varying
%      exposure times, so do that if you want to extract these
%      parameters.  Minimally, you need 2 different exposure times.
%      (Leaving the light off is also a good substitute for changing the
%      exposure time.)
%   3. The gain & noise do not depend on the light level. An important
%      self-consistency check is to take readings at two different light
%      levels and make sure you get the same answers.
%   4. When applying this to PMTs, beware that they are nonlinear except at
%      very high PMT voltage. This then severely restricts the number of
%      photons/pixel before saturation. The best parameters are: high PMT
%      voltage, low laser power, and, if necessary, lots of Kalman
%      averaging.
%   5. In general, you want to avoid having a lot of saturated
%      pixels. While this function can throw out pixels that got
%      saturated, it may be that a saturated pixel also affects its
%      neighbors, which can throw off the computations.
%
% Syntax:
%
% gain = imnoiseanalysis(smm)
%
%   With this syntax, you can extract the gain of the system, the factor
%   that converts quantal units into digital numbers. "smm" is a stackmm
%   object containing your images, or a simple array of images.
%
%
% [gain,bias,noise] = imnoiseanalysis(smm,smmbg)
%
%
%   Also extracts the readout noise and bias of the imaging system. This
%   is possible only if you supply another set of images, smmbg, that
%   contains the "background."  This set of images can be acquired with
%   the light off.
%
%
% [...] = imnoiseanalysis(...,options)
%
%   allows you to control the calculation or get some extra
%   data. "options" is a structure which may have the following fields:
%     maxpixel (default Inf): determines the value returned when the
%       camera/PMT is saturated. Pixels above this value, on any frame,
%       are excluded from the analysis. With "Inf," no pixels are
%       excluded.
%     rescale: an optional vector, rescale(i) is the amount to multiply
%       the ith frame by before analyzing the variance.  This provides an
%       approximate method for compensating for photobleaching.  By
%       default, frames are left unchanged.
%     regression_frac (default [0.1 0.9]): allows you to include only a subset
%       of the pixels in doing the regression, those that are within a
%       certain range of mean values.  The default [0.1 0.9] causes the
%       regression to use pixels in the middle 80% of the range of total
%       intensities.  Set to [0 1] to use all pixels.
%     plotgain (default true): if true, shows the scatter plot of
%       var(pixel) vs. mean(pixel);
%     numbins (default 100): sets the number of bins used in the
%       variance-to-mean plot.
%
% The outputs are as follows:
%   gain is a scalar yielding the gain g.
%   bias is a matrix of the dimensions of a single image.
%   noise is the variance of the (readout) noise, averaged over all
%     pixels of the image.
%
% See also: ANALYZECAMERA.
  
% Copyright 2004 by Timothy E. Holy
  
% Old options:
%     discardmax (default true): if true, throws out any pixel which
%       reaches the maximum value returned by the camera (saturated
%       pixels), as reported by the imrange parameter;
%     plotbias (default 500): if > 0, randomly selects plotbias pixels and
%       plots their mean pixel value vs. integration time;
%     discardmin (default false): if true, throws out any pixel which
%       reaches the minimum value returned by the camera

  smmbg = [];
  have_bias = false;
  options = struct;
  if (nargin > 1)
    if isstruct(varargin{1})
      options = varargin{1};
    elseif (isnumeric(varargin{1}) || isa(varargin{1},'stackmm'))
      smmbg = varargin{1};
      have_bias = true;
      if (nargin > 2)
        options = varargin{2};
      end
    end
  end
  options = imnaoptions(options);

  % If we have inputs that give us the bias, compute it
  if have_bias
    [imdims,n_frames] = imna_getimdims(smmbg);
    colons = repmat({':'},[1 length(imdims)-1]);
    bg = smmbg(:,:,colons{:});
    bg = bg(:,:,:);   % In case bg is 4-dimensional
    bias = mean(bg,3);
  end
  
  % Get image dimensions, and load the image if necessary
  [imdims,n_frames] = imna_getimdims(smm);
  if (n_frames < 2)
    error('We need at least 2 frames to estimate the variance');
  end
  n_dims = length(imdims);
  npix = prod(imdims);
  colons = repmat({':'},1,n_dims);
  m = smm(colons{:},:);  % if necessary, load the whole thing in
  
  % Because we might rescale, we need to compute mean by hand
  fmean = zeros(imdims);
  for frameIndex = 1:n_frames
    mtmp = m(colons{:},frameIndex);
    if isfield(options,'rescale')
      mtmp = mtmp * options.rescale(frameIndex);
    end
    fmean = fmean + double(mtmp);
  end
  fmean = fmean/n_frames;

  fvar = zeros(imdims);
  % Compute variance with a loop, so it doesn't eat memory
  for frameIndex = 1:n_frames
    mtmp = m(colons{:},frameIndex);
    if isfield(options,'rescale')
      mtmp = mtmp * options.rescale(frameIndex);
    end
    fvar = fvar + (double(mtmp) - double(fmean)).^2;
  end
  fvar = fvar/n_frames;
  
  % Keep track of saturated pixels
  flagbad = any(m >= options.maxpixel,n_dims+1);

  % Now compute the bias-subtracted mean
  if have_bias
    fmean = bsxfun(@minus,fmean,bias);
  end

  % Reshape mean & var to "eliminate" spatial info
  fmean = fmean(:);
  fvar = fvar(:);
  flagbad = flagbad(:);
  if (sum(flagbad) > 0)
     warning([num2str(sum(flagbad)) ' saturated pixels found']);
  end
  
  % Calculate gain and noise
  options.biassubtract = have_bias;
  [g,n] = imna_gain_noise(fmean(~flagbad),fvar(~flagbad),options);
  gain = g;
  if have_bias
    noise = n;
  end
  
function options = imnaoptions(options)
%  if ~isfield(options,'plotbias')
%    options.plotbias = 500;
%  end
  if ~isfield(options,'maxpixel')
    options.maxpixel = Inf;
  end
  if ~isfield(options,'plotgain')
    options.plotgain = 1;
  end
  if ~isfield(options,'regression_frac')
    options.regression_frac = [0.1 0.9];
  end
  if ~isfield(options,'numbins')
    options.numbins = 100;
  end
  % Note rescale is deliberately not handled here
  
  
function [g,n] = imna_gain_noise(mtmp,vtmp,options,gainset)
  % gain = slope of variance vs. mean line
  % noise = y-intercept of variance vs. mean line
  [smtmp,sort_order] = sort(mtmp(:));
  svtmp = vtmp(sort_order);
  intensity_range = smtmp(1) + ...
      options.regression_frac * (smtmp(end)-smtmp(1));
  use_for_regression = smtmp >= intensity_range(1) & ...
      smtmp <= intensity_range(2);
  [g,n] = linregress(smtmp(use_for_regression),svtmp(use_for_regression));
  if options.plotgain
    % To make the figure manageable, do some averaging to reduce the
    % number of points
    %splitpoints = round(linspace(1,length(smtmp),options.numbins));
    split_intensity = linspace(smtmp(1),smtmp(end),options.numbins);
    pye = nan(1,options.numbins-1);
    py = pye;
    px = pye;
    for k = 1:length(split_intensity)-1
        binIndex = smtmp >= split_intensity(k) & smtmp < split_intensity(k+1);
        if sum(binIndex)
            px(k) = mean(smtmp(binIndex));
            py(k) = mean(svtmp(binIndex));
            pye(k) = std(svtmp(binIndex))/sqrt(length(binIndex));
        end
    end
    figure
    %h = errorbar(px,py,pye);
    %set(h,'LineStyle','none');
    %h = line([px; px],[py-pye;py+pye],'Color','b','LineWidth',3);
    h = line([px; px],[py-pye;py+pye],'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',1.5,'LineWidth',0.5);
    if options.biassubtract
      xlabel('Bias-subtracted mean pixel value')
    else
      xlabel('Mean pixel value')
    end
    ylabel('Var pixel value')
    pym = px * g + n;
    line(px,pym,'Color','r','LineWidth',2.5);
    if ~options.biassubtract
      n = NaN;
    end
    title(sprintf('Gain=%g grays/quantum, var(noise)=%g grays^2',g,n));
  end
    
function [imdims,n_frames] = imna_getimdims(smm)
  if isa(smm,'stackmm')
    sz = smm.size;
  else
    sz = size(smm);
  end
  if (length(sz) > 2)
    imdims = sz(1:end-1);
    n_frames = sz(end);
  else
    imdims = sz;
    n_frames = 1;
  end
  return;
  
%
% Below is some old code that used CCD integration time as a "better" way
% of estimating bias, so that the gain also applied to the dark
% current. But this is hard to do with PMTs.  So it's abandoned

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
    
  