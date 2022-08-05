function goodchannelso = valchanrun(cycs,vlvusage)
% VALCHANRUN: run validatechannels on a set of ephys structures
%
% This gives the option of selecting time periods with stimulus, and
% examining only those snippets occurring during valve open periods.
% This way, one is likely to be looking at the snippets that are actually
% important for some conclusion.
%
% goodchannels = valchanrun(cycs,vlvusage)
% where
%   cycs is a structure array of type EPHYS.
%   vlvusage indicates how to use stimulus information. If 'ignore', then
%     all snippets are used. If 'all', then all (nonzero) valve numbers are used.
%     Alternatively, a vector of valve numbers can be supplied, and then
%     only periods with those valve numbers are used. Default: 'ignore'.
% and
%   goodchannels (optional) is an ephys structure array of length(cycs), giving
%     the list of "approved" channel numbers in the channels field. If empty,
%     the user is prompted to save this to a .mat file.
%
% See also: VALIDATECHANNELS, EPHYS.

% Determine periods of stimulation
if nargin < 2
  vlvusage = 'ignore';
end
ignoring = 0;
if (ischar(vlvusage) & strcmp(vlvusage,'ignore'))
  ignoring = 1;
else
  goodvlvs = vlvusage;
end
nfiles = prod(size(cycs));
if ~ignoring
  loadindx = ephystofetch(cycs,'stimulus');
  if ~isempty(loadindx)
    cycs(loadindx) = ephysfetch(cycs(loadindx),'stimulus');
  end
  [intervals,identity] = intervalsfromstim(cycs,'open',goodvlvs);
  % Avoid stimulus artifacts
  for i = 1:nfiles
    intervals{i}(:,1) = intervals{i}(:,1) + cycs(i).scanrate/10;
  end
  for i = 1:nfiles
          cycsptmp = ephyssubrange(cycs(i),intervals{i});
          [cycsptmp.tag] = deal(cycs(i).valvelabels{identity{i}});
          cycsp{i} = cycsptmp;
  end
else
  for i = 1:nfiles
    cycsp{i} = cycs(i);
  end
end
% Get snippets
for i = 1:nfiles
    cycsp{i} = ephysfetch(cycsp{i},{'sniptimes','snippets'});
end
% Set up output structure and call validatechannels
[goodchannels(1:nfiles).basefilename] = deal(cycs.basefilename);
for i = 1:length(cycsp)
  tmp = validatechannels(cycsp{i});
  goodchannels(i).channels = tmp;
end
goodchannels = reshape(goodchannels,size(cycs));
% Save, if no output argument
if nargout == 0
  [filename,pathname] = uiputfile('goodchan.mat','Save approved channel numbers');
  if (filename ~= 0)
    save([pathname,filename],'goodchannels');
  end
else
  goodchannelso = goodchannels;
end
