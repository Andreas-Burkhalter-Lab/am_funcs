function [varargout] = spikes_vs_stim(varargin)
% Sorts spikes according to time relative to the stimulus.
%  For each .ssnp file, requires a 'snipstruct.mat' file a in directory with the
%  same name as the .ssnp file, as creaetd by sorthead_from_raw_all. 
%  This function adds to saved snipstruct variable stimulus timing information
%  and sorts snips into pre-stim and post-stim. 
%  Enter 'getstimwave' as an argument to also save the waveform on the
%  specified stimulus channel as a field in snipsstruct.
%  Currently only designed to handle one stimulus time. 
 
% In 15034, stim reaches -8,8 and stimchan noise reaches -1,1
% In 15033, stim reaches +-1.5 and stimchan noise reaches +-0.3
stimthresh = 1;  % maybe make this an optional input argument
stimchan = 63;  % maybe make this an optional input argument

 filenames = dir('*.ssnp');
[placeholder filenames] = cellfun(@fileparts,extractfield(filenames,'name'),'Uniformoutput',false);

for i = 1:length(filenames)
   filenames{i}
   mer = merecmm(strcat(filenames{i},'.merec')); 
   load(strcat(filenames{i},'/snipstruct.mat'));
   clear snipstruct.stimtime snipstruct.tpre snipstruct.tpost snipstruct.pre snipstruct.post
   snipstruct.stimtime = get_stim_time(mer,stimchan,stimthresh);
   snipstruct.stimchan = stimchan;
   snipstruct.stimthresh = stimthresh;
   if ~isempty(varargin) && strcmp(varargin,'getstimwave')
       snipstruct.stimwave = mer([stimchan],[1:mer.nscans]);
   end
   for j = 1:length(snipstruct.channels)
       if ~isempty(snipstruct.sniptimes{j})
           snipstruct.tpre{j} = snipstruct.sniptimes{j}(snipstruct.sniptimes{j} <= snipstruct.stimtime);
           snipstruct.tpost{j} = snipstruct.sniptimes{j}(snipstruct.sniptimes{j} > snipstruct.stimtime);
           snipstruct.pre{j} = snipstruct.snips{j}(:,snipstruct.sniptimes{j} <= snipstruct.stimtime);
           snipstruct.post{j} = snipstruct.snips{j}(:,snipstruct.sniptimes{j} > snipstruct.stimtime);
       end
   end
   save(strcat(filenames{i},'/snipstruct.mat'),'snipstruct')
end

end