function batch_sound2sng_onedir
% BATCH_SOUND2SNG_ONEDIR: processing for large numbers of singing trials
% This function loops over all the .bin files in a directory and creates
% the .sng files. You call it with no arguments.
  binfiles = dirbytime('*.bin');
  nfiles = length(binfiles);
  for idx=1:nfiles
    [pathstr,basename] = fileparts(binfiles{idx});
    if ~exist([basename '.sng'],'file')
      sound2sng(binfiles{idx}, ...
		struct('threshold',0.28,'nfreq',256, ...
		       'freqrange',[25000 110000],...
		       'plot',0),...
		[basename,'.sng'])
    end
  end
