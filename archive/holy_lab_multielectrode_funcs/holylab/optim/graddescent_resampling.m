
  function [p,fval] = graddescent_resampling(func,p0,data,options)
% GRADDESCENT_RESAMPLING: a stochastic gradient descent minimization algorithm
%
% In high-dimensional optimization problems, convergence to a "good"
% minimum is difficult.  Stochastic optimization uses noise to escape
% small local minima.  This function performs stochastic optimization in
% the following way: on each iteration, it chooses a random subset of the
% data and does one round of line minimization.
%
% Syntax:
%   p = graddescent_resampling(func,p0,data,options)
% where
%   func is the function handle to your function to be minimized (which
%     has a syntax [val,grad] = func(p,data), where val is the value and grad
%     is the gradient with respect to the parameters);
%   p0 is your initial guess for the parameters
%   data is a (multidimensional) array of data values
%   options is a structure which may have the following fields:
%     iter_max (default 1000): the maximum number of parameter tries;
%     subset_dimension (default 1): the dimension of data for which you
%       will be choosing the random subset
%     subset_length (default length(p0)): the size of each random subset
%     mu (default 1): the initial factor multiplying the gradient in
%       determining the step size;
%     mu_decfactor (default 5): the factor by which to decrease mu on
%       failed steps;
%     mu_incfactor (default 2): the factor by which to increase mu on
%       successful steps;
%
% See also: FMINUNC, GRADDESCENT.

% Copyright 2009 by Timothy E. Holy

if (nargin < 4)
    options = struct;
end
options = default(options,'iter_max',100000,'subset_dimension',2,'subset_length',500,'mu',1,'mu_decfactor',5,'mu_incfactor',2,'Display',true);

n_data = size(data,options.subset_dimension);

p = p0;
mu = options.mu;
iter = 0;
isdone = false;
fval = []; % history of error
mu_val = []; % history of mu
if options.Display
    figure;
    hax = gca;
    hax2 = SplitVert([0.45 0.55],[1 0 1]);
    axes(hax2(1));
    xlabel('Iteration');
    ylabel('Function value');
    axes(hax2(2));
    xlabel('Iteration');
    ylabel('mu');
end

index = randperm(n_data);
index = index(1:options.subset_length);
cdata = data(index);
[v,g] = func(p,cdata);
fval(1) = v;
while ~isdone
    ptest = p - mu*g;
    g2 = g(:)'*g(:);
    vNew = func(ptest,cdata);
    if (vNew >= v)
        % Bad step, diminish mu
        mu = mu/options.mu_decfactor;
    else
        % Good step, keep it
        p = ptest;
        fval(end+1) = vNew;
        mu_val(end+1) = mu;
        if options.Display
            axes(hax2(1))
            cla;
            line(1:length(fval),fval);

            axes(hax2(2));
            cla;
            line(1:length(mu_val),mu_val);
            drawnow
            end
        % Choose a new random subset of the data
        index = randperm(n_data);
        index = index(1:options.subset_length);
        cdata = data(index);
        [v,g] = func(p,cdata);
        % Try an even bigger mu next time
        mu = mu*options.mu_incfactor;
    end
    iter = iter+1;
    if (iter > options.iter_max)
        isdone = true;
    end
end
fval(end+1) = v;