function err = optmatch(w,K,ctensor,itarget)
% OPTMATCH: designer input currents to mitral cells
%
% This function would generally be supplied to a minimization routine
% such as fminsearch.
%
% Usage:
%   err = optmatch(w,K,ctens,itarget)
% where
%   w is a vector of length ncells giving synaptic weights
%   K is a ncells-by-2 matrix of association constants
%   ctensor is a n-by-m-by-2 tensor, where c(i,j,1) is the concentration
%     of compound #1 at vertex (i,j), and c(i,j,2) is the concentration
%     of compound #2 at this vertex.  For example, if you want to
%     decompose the mitral cell current as a product of f(c1/c2) *
%     g(c1*c2), then ctensor(i,:,:) should correspond to lines of
%     constant c1/c2, and ctensor(:,j,:) should correspond to lines of
%     constant c1*c2.
%     If the sum of concentrations at a vertex is larger than 1, this
%     vertex is not used for generating the best-fit to itarget. 
%   itarget is the target matrix of input currents to the mitral cell.
% and
%   err is the mean-square error 
%
% See also: CTENSRATIO.
  
% Tim Holy, 01-21-2002
[n,m,p] = size(ctensor);
if (p ~= 2)
  error('The concentration tensor must be n-by-m-by-2');
end
err = 0;
for i = 1:n
  for j = 1:m
    % Check to see if the concentration is valid
    if (sum(ctensor(i,j,:)) <= 1)
      err = err + (calccurrent(K,w,ctensor(i,j,:)) - itarget(i,j))^2;
    end
  end
end  
