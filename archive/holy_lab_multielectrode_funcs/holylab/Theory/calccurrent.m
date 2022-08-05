function iinput = calccurrent(K,w,c)
% CALCCURRENT: calculate the input current to a mitral cell
% This function calculates the input current to a "mitral cell",
% given a model of the following form:
%     I = sum_i w_i r_i
% where w_i is a set of synaptic weights and r_i is a set of firing
% rates, which satisfy the following rule:
%     r_i = (sum_j c_j/K_ij)/(1 + sum_j c_j/K_ij).
% Here c_j is the concentration of the jth compound, and K_ij is the
% association constant of the ith neuron for the jth compound.
%
% Syntax:
%   I = calccurrent(K,w,c)
% where
%   K is a numcells-by-numcompounds matrix
%   w is a vector of length numcells
%   c is a vector of length numcompounds
  
c = c(:);
w = w(:);
[ncells,ns] = size(K);
if (length(w) ~= ncells)
  error(['The number of synaptic weights must be equal to the number of' ...
	 ' cells']);
end
if (length(c) ~= ns)
  error(['The number of concentrations must be equal to the number of' ...
	 ' compounds']);
end
Kinv = 1./K;
cnorm = Kinv*c;
r = cnorm./(1+cnorm);
iinput = w'*r;
