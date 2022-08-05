function wav = whis2wav(whis,gap,options)
% whis is a cell array of sparse sonograms during whistles
% gap is a vector of gap lengths between whistles (in # of columns)
% wav is the pitch-shifted reconstructed waveform (normalized) suitable
% for saving to a wav file
  
  % Re-synthesize sparse sonogram
  nwhis = length(whis);
  if (length(gap) + 1 ~= nwhis)
    error('Number of gaps does not match number of whistles');
  end
  nfreq = size(whis{1},1);
  sngpieces = cell(1,2*nwhis-1);
  for i = 1:nwhis-1
    sngpieces{2*i-1} = whis{i};
    sngpieces{2*i} = spalloc(nfreq,gap(i),0);
  end
  sngpieces{end} = whis{end};
  sng = cat(2,sngpieces{:});
  
  % Re-synthesize waveform
  x = sng2sound(sng);
  
  % Pitch-shift waveform
  wav = phasevocoder(x,2*(nfreq-1),16,eps,[0 1]);
  wav = wav/max(abs(wav))/(1+eps);
  
  