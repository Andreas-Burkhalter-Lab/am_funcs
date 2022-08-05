function analysis = seea(varargin) 
% SEAA is a shorthand front for SINGLE_ELECTRODE_EPHYS_ANALYZE
% Usage: analysis = seea( ) -- analyzes current directory only
%        analysis = seea(...,'directory', directory) -- analyzes directory indicated
%        analysis = seea(...,'recursive', ['on'/'off']) -- recurses through directories
%        analysis = seea(..., 'manual_valvelabel', _filename_)
%                   -- will use the input filename to assign valves.  This
%                      is important on recordings where errors in automatic 
%                      valve assignment are evident "filename" must refer
%                      to a function which returns a cell array of strings
%                      including the desired valve labels.
%
% See also SINGLE_ELECTRODE_EPHYS_ANALYZE

analysis = single_electrode_ephys_analyze(varargin);