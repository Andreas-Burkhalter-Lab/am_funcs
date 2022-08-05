function analysis = meea(varargin)
% MEAA is a shorthand for MULTI_ELECTRODE_EPHYS_ANALYZE
% Usage: analysis = meea( ) -- analyzes current directory with default values
%        analysis = meea('param', value, ...) -- sets parameters to paired values
%        analysis = meea(struct('param',value,...)) -- sets parameters to values in input struct fields
%        **NOTE: the "call_func" parameter sets the function (and subsequent
%                calculations) to be performed on your multielectrode data set)
% Output: analysis is a struct containing various subfields calculated by
%         designated call functions
%
% Examples: 
%        a = meea('directory', directory, 'call_func', @meea_deltar)
%            -- analyzes directory indicated using meea_deltar.m to calculate delta_r
%        a = meea(..., 'manual_valvelabel', _filename_)
%                   -- will use the input filename to assign valves.  This
%                      is important on recordings where errors in automatic 
%                      valve assignment are evident "filename" must refer
%                      to a function which returns a cell array of strings
%                      including the desired valve labels.
%
% See also MULTI_ELECTRODE_EPHYS_ANALYZE

analysis = multi_electrode_ephys_analyze(varargin);