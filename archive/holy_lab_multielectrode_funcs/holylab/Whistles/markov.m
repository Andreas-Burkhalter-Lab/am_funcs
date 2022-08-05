function [ustates,p,n] = markov(seq)
% MARKOV: model a sequence using a Markov model
% Syntax:
%   [ustates,ptrans] = markov(seq)
%   [ustates,ptrans,ntrans] = markov(seq)
% where
%   seq is a vector containing the states of the model (e.g.,
%     0, 1, 2, ..., or 'a', 'none', ...)
% and
%   ustates is a vector containing the unique states;
%   ptrans is an nstates-by-nstates matrix, ptrans(i,j) is the observed
%     frequency of transitions from ustates(i) to ustates(j);
%   ntrans is an nstates-by-nstates matrix, ntrans(i,j) is the number of
%     observed transitions from ustates(i) to ustates(j) (in other words,
%     ptrans is just a normalized version of ntrans).

% Copyright 2005 by Timothy E. Holy
  
  [ustates,tmp,stateindx] = unique(seq);
  nstates = length(ustates);
  n = zeros(nstates);
  for i = 1:length(stateindx)-1
    n(stateindx(i),stateindx(i+1)) = n(stateindx(i),stateindx(i+1)) + 1;
  end
  ntot = sum(n,2);
  p = n ./ repmat(ntot,1,nstates);
  