%  register_barrierpenalty_composition
%  
%  The actual calculation is implemented in register_barrierpenalty_composition.cpp
%  which calculates the edge barrier penalty for register_block_penalty.m.
% 
%  The usage of this mex code is
%  [BPsum, BPGrad] = register_barrierpenalty_composition(unew(gridIndex,:),mismatch{gridIndex});
%  where
%  BPsum is the total edge barrier penalty,
%  BPGrad is the gradient of BPsum in this mismatch cell.

%  Copyright 2012 Jian Wang




