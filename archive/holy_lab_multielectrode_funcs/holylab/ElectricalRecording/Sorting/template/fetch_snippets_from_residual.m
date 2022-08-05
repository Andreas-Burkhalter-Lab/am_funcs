function snips = fetch_snippets_from_residual(amplitudes,templates,times,memmres,channels,sniprange)
% Syntaxes:
%   snips = fetch_snippets_from_residual(amplitudes,templates,times,memmres,channels,sniprange)
% templates must be a templateLen-by-n_templates matrix
  snips = templates * amplitudes;
  % Now add in the residual
  if (nargin > 2)
    snipLen = diff(sniprange)+1;
    nSnips = length(times);
    snipRes = fetch_snippets_from_merecmm(memmres,channels,times,sniprange);
    % We have to change it to time,channel,snips format & then unweave the
    % first two indices
    snipRes = permute(snipRes,[2 1 3]);
    snipRes = reshape(snipRes,[length(channels)*snipLen nSnips]);
    snips = snips + snipRes;
  end
  