function ephys = ephyssimulatedspikes(eparams,rparams)
% EPHYSSIMULATEDSPIKES: generate "responses" in ephys structure
% Syntax:
%   ephys = ephyssimulatedspikes(eparams,rparams)
% where
%   eparams is a structure with the following fields:
%     nrepeats is the number of repeats;
%     stimulus is a 2-by-n matrix, stimulus(1,:) is the sequence of valve
%       states and stimulus(2,:) contains the transition times, measured in
%       seconds (note stimulus(2,end)-stimulus(2,1) specifies the duration of
%       the window);
%     scanrate is the fake scan rate (optional, default 10000);
%     tag is a string indentifying this "stimulus";
%   rparams is a structure with the following fields:
%     ratefunc is a function handle or string with the name of the rate
%       function that will be used to generate the fake spikes, with the syntax
%           r = ratefunc(t,stimulus_secs,rparams)
%       where t is a vector of times, stimulus_secs is a stimulus matrix
%       where the times have been converted into seconds, and r can be a
%       ncells-by-ntimepoints matrix of rates if multiple cells are to be
%       simulated;
%     + any parameters that need to be passed to the rate function to
%       specify the firing rate(s);
%
% and
%   ephys is the output ephys structure.
%
% See also: RATE1STORDER.

% Set up the appropriate output structure
if ~isfield(eparams,'scanrate')
  eparams.scanrate = 10000;
end
stimulus_secs = eparams.stimulus;
ephys.stimulus = eparams.stimulus;
ephys.stimulus(2,:) = ephys.stimulus(2,:)*eparams.scanrate;
ephys.tag = eparams.tag;
ephys.scanrate = eparams.scanrate;
ephys.scanrange = eparams.stimulus(2,[1 end])*eparams.scanrate;
ephys.toffset = -diff(eparams.stimulus(2,[1 2]));  % First transition time

% Get the firing rates & finish preparing structure
ntimepoints = 10000;
trange = stimulus_secs(2,[1 end]);
t = linspace(trange(1),trange(2),ntimepoints);
dt = diff(trange)/(ntimepoints-1);
r = feval(rparams(1).ratefunc,t,stimulus_secs,rparams);
ncells = size(r,1);
ephys.channels = 1:ncells;
ephys = repmat(ephys,1,eparams.nrepeats);

% Compute \int dt r(t); create an "inverse map," so we can place spikes
% uniformly in an interval and then map them to real times.
intr = dt*cumsum(r,2);
for j = 1:ncells
  Nspikes = intr(j,end);
  [s,sindx] = sort([linspace(0,Nspikes,ntimepoints) intr(j,:)]);
  tmap(j,:) = find(sindx <= ntimepoints) - (0:ntimepoints-1);
end
tmap = (tmap/ntimepoints)*diff(trange) + trange(1);

% Create simulated spikes, using Poisson statistics.
for i = 1:eparams.nrepeats
  for j = 1:ncells
    Nspikes = intr(j,end);  % This is the _expected_ number
    twait = expdev(3*ceil(Nspikes)+40);  % interspike intervals
    tspike = cumsum(twait);
    if (max(tspike) < Nspikes)
      error('Didn''t generate enough spikes');
    end
    tspike = tspike(find(tspike < Nspikes));
    % Now map to real times by indexing into tmap
    tspike = tspike/Nspikes * (ntimepoints-1);
    tspfl = ceil(tspike);
    tspfr = tspike-tspfl;
    treal = (1-tspfr).*tmap(j,tspfl) + tspfr.*tmap(j,tspfl+1); % interpolate
    ephys(i).sniptimes{j} = treal*eparams.scanrate;
  end
end

function t = expdev(N)
x = rand(1,N);
indx0 = find(x == 0);
while ~isempty(indx0)
  x(indx0) = rand(1,length(indx0));
  indx0 = indx0(find(x(indx0) == 0));
end
t = -log(x);
