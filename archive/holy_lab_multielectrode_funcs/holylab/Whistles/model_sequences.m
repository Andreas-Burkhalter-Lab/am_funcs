function out = model_sequences(class,t,options)
% MODEL_SEQUENCES: compare distribution of sequences to Markov models
% Syntax:
%   s = model_sequences(class,t,options)
% where
%   class is the syllable type;
%   t is the time at which each chirp started (may be left empty);
%   options (optional) is a structure with the following fields:
%     seqlength (default 5): the length of sequences to consider
% and
%   s is an output structure with the following fields:
%     word: a list of all "words" (sequences of symbols) of length
%       seqlength;
%     nword: a vector containing the number of observations of each word;
%     nfreq: a prediction for the number of words that would have been
%       observed if the sequence were random, based on the relative
%       frequency of each symbol;
%     nmark: a prediction for the number of words that would have been
%       observed based on a Markov model, i.e., using the probability of
%       the first symbol and then the transition probabilities for each
%       2-symbol transition;
%     state: a list of the unique symbols;
%     state_n: the number of times each state was observed;
%     state_freq: the probability of each symbol;
%     tprob: a matrix, tprob(i,j) is the probability that a symbol of
%       index i will be followed by a symbol of index j.
%
% Note: a better way to handle gaps is by introducing a gap state to the
% model. Hence, don't use the results of this function for anything
% except word analysis---the "state" outputs are basically legacy. See
% MARKOV (& introduce a gap state).
%
% See also: SEQUENCES, MARKOV.
  
% Copyright 2005 by Timothy E. Holy
  if (nargin < 3)
    options = struct;
  end
  options = ms_parse(options);
  % Calculate the distribution of states
  state = unique(class);
  nstates = length(state);
  state_freq = nan(1,nstates);
  for i = 1:nstates
    state_n(i) = length(find(class == state(i)));
  end
  state_freq = state_n/length(class);
  % Calculate the transition probabilities in Markov model
  ops2 = options;
  ops2.seqlength = 2;
  [seqs2,n2] = sequences(class,t,ops2);
  tprob = nan(nstates,nstates);
  for i = 1:nstates
    indx = find(seqs2(:,1) == state(i));  % Find all transitions from
                                          % given state
    snorm = sum(n2(indx));
    for j = 1:length(indx)
      stateindx = find(state == seqs2(indx(j),2));  % Figure out dest. state
      tprob(i,stateindx) = n2(indx(j))/snorm;
    end
  end
  % Measure distribution of "words"
  [word,n] = sequences(class,t,options);
  nuseqs = length(n);
  ntotseqs = sum(n);
  wordprob = n/ntotseqs;
  % Calculate expected distribution of words from likelihood of
  % individual symbols
  freqprob = ones(1,nuseqs);
  for i = 1:nuseqs
    for j = 1:length(word(i,:))
      k = find(state == word(i,j));
      freqprob(i) = freqprob(i)*state_freq(k);
    end
  end
  % Calculate expected distribution of words from Markov transition model
  markprob = ones(1,nuseqs);
  for i = 1:nuseqs
    k = find(state == word(i,1));
    markprob(i) = state_freq(k);   % Probability of first symbol
    for j = 1:length(word(i,:))-1
      k1 = find(state == word(i,j));
      k2 = find(state == word(i,j+1));
      markprob(i) = markprob(i)*tprob(k1,k2);  % Prob of transition
    end
  end
  out.word = word;
  out.nword = n;
  %out.wordprob = wordprob;
  out.nfreq = freqprob*ntotseqs;
  out.nmark = markprob*ntotseqs;
  out.state = state;
  out.state_n = state_n;
  out.state_freq = state_freq;
  out.tprob = tprob;
  
  
function ops = ms_parse(ops)
  if ~isfield(ops,'seqlength')
    ops.seqlength = 5;
  end
  