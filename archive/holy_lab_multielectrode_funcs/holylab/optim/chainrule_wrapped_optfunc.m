function [varargout] = chainrule_wrapped_optfunc(p,optfunc,pin2pout,gpout2gpin)
% CHAINRULE_WRAPPED_OPTFUNC: a default wrapper for chain rule optimization
%
% For optimization problems involving large numbers of variables, it
% is sometimes desirable to be able to test a solution satisfying a more
% restricted form.  For example, in a problem which is naturally expressed
% in terms of a "field" of variables (i.e., a numeric array), one might be
% interested in a model in which the field is restricted to be in the form
% of a particular parametric model. Rather than implementing two versions
% of your optimization function, you can implement the general "field"
% version and then use the chain rule to convert it into a model over
% particular parametric forms.  You just need two function handles, one
% that computes the field from the parameters, and a second that convert
% gradients with respect to the field into gradients with respect to the
% parameters.
%
% General usage pattern: imagine that the "full" problem operates on a
% multidimensional array A whose values are to be optimized.  Therefore, it
% is assumed that you have written a function with the following syntaxes
% available:
%   val = optfunc(A);
%   [val,grad] = optfunc(A);
% Now suppose you want to restrict A to have a certain parametric form,
% here assumed for the purposes of illustration to be Gaussian. You
% want to optimize the ampltidue, mean, and covariance of the
% Gaussian. You will need to define two function handles like this:
%   [pin2pout,gpout2gpin] = chainrule_gaussian(x,insertionFlag)
% (see CHAINRULE_GAUSSIAN for details)
% and then create the final function like this:
%   goptfunc = @(p) chainrule_wrapped_optfunc(p,optfunc,pin2pout,gpout2gpin);
% Then you can use goptfunc in the minimization program, e.g.,
%   p = fminunc(goptfunc,p0,optimset('GradObj','on'))
% where p and p0 are vectors specifying the parameters of the gaussian
% (see below).
%
% Details on syntax:
%   [val,grad] = chainrule_wrapped_optfunc(p,optfunc,pin2pout,gpout2gpin)
% where
%   p is a vector of "input parameters"
%   optfunc is an optimization function that works on "output parameters"
%   pin2pout converts "input parameters" (pin) into "output parameters"
%     (pout)
%   gpout2gpin converts gradients with respect to the "output parameters"
%     into gradients with respect to the "input parameters."
% and
%   val and grad are the function value and gradient, respectively, with
%     respect to the "input parameters" p.
%
% A fully-worked example is described in CHAINRULE_GAUSSIAN.
%
% See also: CHAINRULE_GAUSSIAN.

% Copyright 2009 by Timothy E. Holy

% Note that this function exists more for the "help text" than in doing
% anything terribly complicated!

  pout = pin2pout(p);
  varargout = cell(1,nargout);
  [varargout{:}] = optfunc(pout);
  if (nargout > 1)
    varargout{2} = gpout2gpin(varargout{2});
  end
end