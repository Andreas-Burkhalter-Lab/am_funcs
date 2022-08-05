function perf = spimopt(x,opt_params,fieldnames)
% SPIMOPT: SPIM performance optimization criterion
% Syntax:
%   perf = spimopt(x,opt_params,fieldnames)
% where
%   x is a vector of parameters being modified;
%   opt_params are a set of optical parameters accepted by
%     spimtraceconfig;
%   fieldnames is a cell array of strings corresponding to the names of
%     fields in opt_params; their order corresponds to the order of
%     values in x and replaces the corresponding value in opt_params
%     before calling spimtraceconfig;
% and
%   perf is a scalar output corresponding to the "performance" that you
%     want to optimize, e.g., the width of the waist at the minimum, or
%     perhaps some integrated width over the fieldsize.
%
% See also: SPIM_PCX, SPIMTRACECONFIG.
  
  for i = 1:length(fieldnames)
    opt_params.(fieldnames{i}) = x(i);
  end
  sigma = spimtraceconfig(opt_params);
  perf = min(sigma);
  