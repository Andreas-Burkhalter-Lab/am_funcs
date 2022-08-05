function merec2vlv(outfile,merecfiles,options)
% MEREC2VLV: create a vlv file from merec files
% Syntax:
%   merec2vlv(outfile,merecfiles)
%   merec2vlv(outfile,merecfiles,options)
%
% where
%   outfile:    a string containing the name of the output vlv file
%   merecfiles: a cell array of merec file names (or a string if you're
%               processing a single file). 
%   options:  an options structure with the following fields:
%      nocheck: if true, reads stimulus information directly from the header,
%               skipping the step of checking the analog stimulus copy in the raw
%               waveform. (You can use this option to generate a vlv file from
%               snippet or envelope files.)
%      feedback: the feedback channel (@see timestim())
%      timestim_calc_vs_intended_tolerance: see timestim.m
%   (see timestim.m for more)
% 
% @history:
%    9/29/2003: use parse_stim_seq() to get stimmatrix
%    8/10/2006: added handling for new timestim option (RCH)
%   10/27/2006:  " (RCH)
%
% See also: autosnip2, timestim
  
  if ~iscell(merecfiles)
    if ischar(merecfiles)
      merecfiles = {merecfiles};
    else
      error(['The second argument must be a string or a cell array of' ...
             ' strings']);
    end
  end
  nfiles = length(merecfiles);
  if (nargin < 3)
    options = struct;
  end
  
  for i = 1:nfiles
    h = readheader(merecfiles{i});
    
    if(is_robot_stim(h))
       options.nocheck=0; % force check for robot deliveried stimuli
       nostim=0; % has stimuli
    else
       % stimmatrix = eval(['[' key2value(h.wholeheader, 'stimulus sequence') ']']);
       stimmatrix = parse_stim_seq(h);
       nostim = 0;
       if isempty(stimmatrix)
          nostim = 1;
          stimtemp = [0 0; 0 h.nscans];
       else
          stimtemp = [stimmatrix(:,1) round(stimmatrix(:,2)*h.scanrate)];
       end
    end
    
    if (~isfield(options,'nocheck') || options.nocheck == 0)
      if ~nostim
        % Load raw waveform data and find the exact transition times,
        % storing the corrected value in stimtemp
        timestimoptions = struct;
        if isfield(options,'feedback')
            timestimoptions.feedback = options.feedback;
        end
        if isfield(options,'timestim_calc_vs_intended_tolerance') 
            timestimoptions.calc_vs_intended_tolerance = options.timestim_calc_vs_intended_tolerance;
        end
        if isfield(options,'timestim_wnds')
            timestimoptions.wndL = options.timestim_wnds(1);
            timestimoptions.wndR = options.timestim_wnds(2);
        end
        if isfield(options,'no_intended_times')
            timestimoptions.no_intended_times = options.no_intended_times;
        end
        if isfield(options,'total_override')
            timestimoptions.total_override = options.total_override;
        end
           stimtemp=timestim(merecfiles{i}, timestimoptions);
        %stimtemp=[stimmatrix(:,1) tRealTime*h.scanrate];
      end
    end
    stim{i} = stimtemp';
    [pathtemp,basetemp] = fileparts(merecfiles{i});
    basename{i} = basetemp;
  end
  
  writevlv(outfile,basename,stim);

