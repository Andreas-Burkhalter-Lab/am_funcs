function whiscorr(p,options)
% WHISCORR: compute temporal correlations across whistles
%
% This function computes autocorrelations in the sequence of individual
% whistles for a number of parameters which describe each whistle.  These
% autocorrelations are plotted.
%
% Syntax:
%   whiscorr(p,options)
% where
%   p is a parameters structure, computed from WHISPARAMS;
%   options is an optional structure controlling the presentation of
%     results, using the following fields:
%       shuffle: if true, also compute correlations when parameters are
%         shuffled (TODO: value determines the number of shuffles);
%       fieldnames: a cell array of parameter names (parameters from
%         WHISPARAMS) determining which parameters are examined;
%       logflag: a vector of length(fieldnames), if true the log of the
%         corresponding field is taken before computing correlations;
%       mode: correlations between whistles are examined as a function of
%         their timing (if mode = 'time') or simply in terms of their
%         order (mode = 'number'). Default 'time'.
%       corlen: the length of correlations.  In 'time' mode, this is in
%         seconds, while in 'number' mode it simply counts the number of
%         intervening whistles (default 4);
%       nbins: used only in 'time' mode, sets the number of bins used in
%         averaging the correlation (default 50);
%       meanonly: (time mode) if true, does not show a scatter plot of the
%         individual products
%
% See also: WHISPARAMS.

% Copyright 2004 by Timothy E. Holy
  
% Set up defaults
if (nargin < 2)
  options = struct;
end
if isempty(p.wt)
  return;
end
%fieldnames = {'dt','pow','meanf','meanv','varv','ratv','meana','vara','rata'};
fields = {'maxjump','dt','pow','dtheta','haslj'};
if isfield(options,'fieldnames')
  fieldnamesall = options.fieldnames;
end
exceptionfields = intersect(fields,{'peakfreq'});
scalarfields = setdiff(fields,exceptionfields);
fields = {scalarfields{:},exceptionfields{:}}; % Re-order bcz alphabetization
nfields = length(fields);
logflag = zeros(1,nfields);
if isfield(options,'logflag')
  logflag = options.logflag;
end
if (length(logflag) < nfields)
  logflag(end+1:nfields) = 0;
end
if (length(logflag) ~= nfields)
  error('logflag and fieldnames don''t match');
end
for i = 1:nfields
  if logflag(i)
    labelname{i} = ['log(' fields{i} ')'];
  else
    labelname{i} = fields{i};
  end
end
%labelname{end+1} = 'fvst';
corlen = 4;
if isfield(options,'corlen')
  corlen = options.corlen;
end
mode = 'time';
if isfield(options,'mode')
  mode = options.mode;
end
if ~strcmp(mode,'number')
  nbins = 50;
  if isfield(options,'nbins')
    nbins = options.nbins;
  end
end
shuffle = 0;
if (isfield(options,'shuffle'))
  shuffle = options.shuffle;
end


% Copy over data from scalar fields
for i = 1:length(scalarfields)
  if logflag(i)
    v(i,:) = log(p.(scalarfields{i})(p.indx));
  else
    v(i,:) = p.(scalarfields{i})(p.indx);
  end
end
if shuffle
  shufindx = randperm(size(v,2));
end

% Calculate autocorrelations
if strcmp(mode,'number')
  % Correlating by sequence
  for i = 1:length(scalarfields)
    actmp = xcorr(v(i,:),corlen,'unbiased');
    actmp = actmp(corlen+1:end);
    ac(i,1:length(actmp)) = actmp/actmp(1);
  end
else
  % Correlating by time separation
  tbins = (1:nbins)*(corlen/nbins);
  [binmean,binsem,npb,tacindx,binlist,tac,ac] = autocorrpulse(p.wt(1,p.indx),tbins,v);
  if shuffle
    [shufmean,shufsem] = autocorrpulse(p.wt(1,p.indx),tbins,v(:,shufindx));
  end
end

% Now compute correlations in the "exception fields"
% PeakFreq. vs. time
if ~isempty(strmatch('peakfreq',exceptionfields,'exact'))
  if strcmp(mode,'number')
    warning('peakfreq correlation not implemented for number mode')
  else
    v = p.peakfreq(p.indx);
    [binmean(end+1,:),binsem(end+1,:),ac(end+1,:)] = wavecorr(v,tacindx,binlist,npb);
    if shuffle
      [shufmean(end+1,:),shufsem(end+1,:)] = wavecorr(v(shufindx),tacindx,binlist,npb);
    end
  end
end
if (nargout == 0)
  if strcmp(mode,'number')
    hlines = plot(ac');
    legend(hlines,labelname{:})
  else
    dims = CalcSubplotDims(nfields+1);
    if ~isfield(options,'meanonly')
      for i = 1:nfields
        subplot(dims(1),dims(2),i);
        scatter(tac,ac(i,:),'.');
        line(tbins,binmean(i,:),'Color','r');
        title(labelname{i})
      end
      figure;
    end
    %normbinmean = binmean ./ repmat(binmean(:,1),1,size(binmean,2));
    normbinmean = binmean ./ repmat(max(binmean,[],2),1,size(binmean,2));
    %     hlines = plot(normbinmean');
    %     legend(hlines,labelname{:});
    for i = 1:nfields
      subplot(dims(1),dims(2),i);
      errorbar(tbins,binmean(i,:),binsem(i,:),'r');
      title(labelname{i})
      if shuffle
        hold on
        errorbar(tbins,shufmean(i,:),shufsem(i,:),'k--');
        hold off
      end
    end
  end
else
  
end

function [mn,sem,ac] = wavecorr(v,tacindx,binlist,npb)
  nac = size(tacindx,2);
  ac = zeros(1,nac);
  for j = 1:nac
    v1 = v{tacindx(1,j)}; v2 = v{tacindx(2,j)};
    nmin = min(length(v1),length(v2));
    nmax = max(length(v1),length(v2));
    ac(j) = sum((v1(1:nmin)-mean(v1(1:nmin))) ...
                .*(v2(1:nmin)-mean(v2(1:nmin))));% ...
               %*(nmin/nmax); % Equiv. to padding w/ zeros (?)
    %ac(j) = sum(v1(1:nmin).*(v2(1:nmin)));
    %ac(j) = range(v1(1:nmin))*range(v2(1:nmin));
  end
  % Average over time bins
  mn = nan(1,length(binlist));
  sem = mn;
  for j = 1:length(binlist)
    if ~isempty(binlist{j})
      mn(j) = mean(ac(binlist{j}));
      sem(j) = std(ac(binlist{j}))/sqrt(npb(j));
    end
  end
  