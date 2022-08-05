function p = whisparams(whis,header,wt,options)
% WHISPARAMS: compute parameters describing vocalizations
% The main tool is to replace the sonogram by it's dominant frequency at
% each time.  Several quantities are then used to characterize this function.
% Syntax:
%   p = whisparams(whis,header,wt,options)
% where
%   whis is a cell array of whistle snippets (from sparse sonogram)
%   header is the .sng file header
%   wt is the 2-by-n matrix of whistle start & stop times
%   options primarily controls the method used to estimate the dominant
%     frequency; it takes the following fields:
%       freqrange: the range of frequencies (in Hz) retained in all
%         computations (out-of-range frequencies discarded);
%       medfilter: if true, the dominant frequency is
%         median-filtered. This parameter then sets the size of the
%         median filter (minimum size 3).
%       discardends: if true, the beginning and end points of the
%         freq. vs. time curve are discarded; the value sets the number
%         of points discarded at each end.
%       boxfilt: if set, the freq. vs. time curve is smoothed with a
%         boxcar filter before computing other measures (value sets the
%         filter size).
%       skipzeros: concatenate across time periods with no power (has no
%         effect on the "real" parameters, which are computed from the peak
%         frequency, for which this happens automatically. However, it
%         does influence the computation of the "meanfrequency" vs. time,
%         which historically was used before switching to peakfrequency.)
%
% The output p contains a fairly large number of fields.  The only
% reasonable documentation for these fields is in the source code.
%
% See also: WHISCORR, WHISTIMES, WHISSNIP, WHISPARAMSSUBSET.
  
% Copyright 2004 Timothy E. Holy.
  
  nwhis = length(whis);
  df = header.scanrate/2/header.nfreq;
  if isfield(options,'freqrange')
    f = linspace(0,header.scanrate/2,header.nfreq);
    findx = find(f >= options.freqrange(1) & f <= options.freqrange(2));
  else
    findx = 1:header.nfreq;
  end
  
  p.wt = wt;
  % Duration:
  p.dt = diff(wt);
  % Gap size:
  p.gap = wt(1,2:end)-wt(2,1:end-1);
  
  % Calculate mean & peak frequencies:
  for i = 1:nwhis
    % Truncate out-of-range frequencies, without changing the size of the
    % matrix
    curwhis = whis{i};
    junkindx = setdiff(1:header.nfreq,findx);
    curwhis(junkindx,:) = 0;

    pow = abs(curwhis).^2;
    totPower = sum(pow);  % Total power at each time point
    p.pow(i) = full(sum(totPower));
    nzindex = find(totPower);
    %Mean
    [x,y,s] = find(curwhis);
    spfreq = sparse(x,y,x,size(curwhis,1),size(curwhis,2)); % Frequency
                                                            % weight matrix
    meanfreq = sum(spfreq.*pow)*df;
    meanfreq(nzindex) = meanfreq(nzindex)./totPower(nzindex);
    % Dominant frequency ( = pitch = "peakfreq")
    [maxp,maxpi] = max(pow); % finds max in each column
    [x,y,s] = find(maxp); % finds non-zero values in maxp, and dumps row and column
    peakfreq = maxpi(y)*df; % 
    pfpow = full(maxp(y));
    spow = full(totPower(y));
    pforig = peakfreq;
    % Be more selective
    if isfield(options,'skipzeros')
      meanfreq = meanfreq(nzindex);
      % Not nec for peakfreq, as this by necessity skips empty columns
    end
    if isfield(options,'medfilter')
      if (options.medfilter < 3)
        options.medfilter = 3;
      end
      if (~mod(options.medfilter,2))
        options.medfilter = options.medfilter+1; % Make it odd
      end
      medfltend = floor((options.medfilter-1)/2);
      % Use medfilt rather than medfilt1 to handle edges properly
      meanfreq = medfilt(meanfreq,options.medfilter);
      peakfreq = medfilt(peakfreq,options.medfilter);
    end
    if isfield(options,'discardends')
      endsz = options.discardends;
      meanfreq = meanfreq(1+endsz:end-endsz);
      peakfreq = peakfreq(1+endsz:end-endsz);
      pfpow = pfpow(1+endsz:end-endsz);
      spow = spow(1+endsz:end-endsz);
      pforig = pforig(1+endsz:end-endsz);
    end
    p.meanfreq{i} = meanfreq;
    p.peakfreq{i} = peakfreq;
    p.pfpow{i} = pfpow;
    p.spow{i} = spow;
    p.pforig{i} = pforig;
  end
  
  % Variable declarations
  % (will calculate shortly)
  p.adjf = cell(1,nwhis);  % "Adjusted" (smoothed) freq. vs. time
  p.dadjf = cell(1,nwhis); % Derivative of the above
  p.maxjump = zeros(1,nwhis)+nan;% Maximum abs. value of derivative
  p.thetarange = zeros(2,nwhis)+nan;  % phase start and end
 
  % We're going to need to filter the peak frequency waveform for some of
  % these analyses
  if isfield(options,'boxfilt')
    b = ones(1,options.boxfilt);
    b = b/sum(b);
    a = zeros(size(b));
    a(1) = 1;
  else
    a = [1 0];
    b = [1 0];
  end
  
  % Calculate the maximum jump (good index for separating steps)
  % and a smoothed representation
  for i = 1:nwhis
    freqvst = p.peakfreq{i};
    if (length(freqvst) > 1)
      p.maxjump(i) = max(abs(diff(freqvst)));
      if (length(freqvst) >= 3*length(b))
        freqvst = filtfilt(b,a,freqvst);
        p.adjf{i} = freqvst;
        p.dadjf{i} = diff(freqvst);
      end
    end
  end
  
  % Compute some overall statistics on the frequencies and slopes (needed
  % for later normalization)
  fall = cat(2,p.adjf{:});
  dfall = cat(2,p.dadjf{:});
  stdf = std(fall);
  stddf = std(dfall);
  
  % Compute phase information
  for i = 1:nwhis
    if ~isempty(p.adjf{i})
      f = p.adjf{i};
      [theta,rho] = cart2pol((f(1:end-1)-mean(f))/stdf,p.dadjf{i}/stddf);
      % Have to transform theta so it "wraps around" 2pi boundaries
      dtheta = diff(theta);
      dtheta = mod(dtheta,2*pi);
      subindx = find(dtheta > pi);
      dtheta(subindx) = dtheta(subindx) - 2*pi;
      theta = [0,cumsum(dtheta)]+theta(1); % Now it's fixed up
      %p.thetarange(:,i) = theta([1 end])';
      p.thetarange(:,i) = [min(theta);max(theta)];
    end
  end
  p.dtheta = diff(p.thetarange);
  p.indx = find(~isnan(p.dtheta));
  
  % Compute jump sizes vs. frequency
  p.dpf = cell(1,nwhis);
  p.pfs = cell(1,nwhis);
  for i = 1:nwhis
    if (length(p.peakfreq{i}) > 1)
      p.dpf{i} = diff(p.peakfreq{i});
      p.pfs{i} = p.peakfreq{i}(1:end-1); % "Starting" peak freq
    end
  end
  %p.dpf = cat(2,dpf{:});
  %p.pfs = cat(2,pfs{:});
  %p.pft = p.dpf + p.pfs;  % "Target" peak freq

  % Classify by jumps
  p.haslj = zeros(1,length(p.pfs));
  p.hashj = p.haslj;
  for i = 1:length(p.pfs)
    [ljdtmp,ljutmp,hjtmp] = whisjclassify(p.pfs{i},p.dpf{i});
    if ~isempty([ljdtmp ljutmp])
      p.haslj(i) = 1;
    end
    if ~isempty(hjtmp)
      p.hashj(i) = 1;
    end
  end
  return;


  p.meanf = nan(1,nwhis); % The mean (peak) frequency, averaged across whistle
  p.varf = nan(1,nwhis);  % The variance of the frequency
  p.v = cell(1,nwhis);    % The "frequency velocity"
  p.meanv = nan(1,nwhis); % The "freq. vel." averaged across the whistle
  p.varv = nan(1,nwhis);  % The variance of the freq. vel.
  p.powv = nan(1,nwhis);  % The "freq. vel. power" (2nd moment of the velocity)
  p.a = cell(1,nwhis);    % The "frequency acceleration"
  p.meana = nan(1,nwhis); % etc.
  p.vara = nan(1,nwhis);
  p.powa = nan(1,nwhis);
  % Many of the vocalizations would appear to be well-modelled as a
  % particle moving under constant force (e.g., trajectories in
  % gravity). For each whistle, let's compute the parameters of the
  % relevant model
  for i = 1:nwhis
    freqvst = p.peakfreq{i};
    if (length(freqvst) >= 3*length(b))
      freqvst = filtfilt(b,a,freqvst);
      p.adjf{i} = freqvst;
      % Value measures
      p.meanf(i) = mean(freqvst);
      p.varf(i) = mean((freqvst - p.meanf(i)).^2);
      % Sweep velocity parameters
      vel = diff(freqvst);
      p.v{i} = vel;
      p.meanv(i) = mean(vel);
      p.varv(i) = mean((vel - p.meanv(i)).^2);
      p.powv(i) = mean(vel.^2);
      % Acceleration parameters
      accel = diff(freqvst,2);
      p.a{i} = accel;
      p.meana(i) = mean(accel);
      p.vara(i) = mean((accel - p.meana(i)).^2);
      p.powa(i) = mean(accel.^2);
    end
  end
  p.indx = find(~isnan(p.vara));
  p.pow = full(p.pow);
  % Normalized velocity and acceleration fluctuations
  p.ratv = nan(size(p.varv));
  p.ratv(p.indx) = p.varv(p.indx)./p.powv(p.indx);
  p.rata = nan(size(p.vara));
  p.rata(p.indx) = p.vara(p.indx)./p.powa(p.indx);
  p.ratav = nan(size(p.vara));
  p.ratav(p.indx) = p.vara(p.indx)./p.powv(p.indx);
  