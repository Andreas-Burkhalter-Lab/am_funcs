function ctensor = ctensratio(ratmin,ratmax,nrat,prodmin,nprod)
% CTENSRATIO: make the ctensor matrix for concentration ratios & products
%
% This function makes the ctensor matrix for functions such as
% PRODCURRENT.  See this function for an explanation of purpose.
%
% Usage:
%   ctensor = ctensratio(ratmin,ratmax,nrat,prodmin,nprod)
% where
%   ratmin is min(log10(c1/c2))
%   ratmax is max(log10(c1/c2))
%   nrat is the number of lines spanning this range (logarithmically
%     spaced)
%   prodmin is min(log10(c1*c2))
%   nprod is the number of lines between here and c1=c2=1
%     (logarithmically spaced).
% Note that some of the combinations of concentrations that are generated
%   by this function are not physically realizable.  However, the
%   "error-checking" is done by PRODCURRENT, as this simplifies
%   representation considerably.
%
% See also: PRODCURRENT.

% Tim Holy, 01-21-2002.
ctensor = zeros(nrat,nprod,2);
crat = logspace(ratmin,ratmax,nrat);
cprod = logspace(prodmin,0,nprod);
ctensor(:,:,1) = sqrt(crat'*cprod);
ctensor(:,:,2) = sqrt((1./crat')*cprod);
