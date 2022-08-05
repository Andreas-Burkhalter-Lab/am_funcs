function filterint16
% FILTERINT16: filter int16 (electrophysiology) data efficiently.
% Y = FILTERINT16(B,A,X,CHANINDEX) filters the data in the rows of X with
% the filters described by A and B to create the filtered data Y. Only
% those rows indexed in CHANINDEX are processed; consequently, Y is
% of size length(CHANINDEX)-by-width(X).
%
% If filters A and B are vectors, then they are applied to each row of X
% indexed by CHANINDEX.  If A and B are matrices, then each _column_ of
% the filters is applied to the corresponding indexed row of X. If A is
% empty, it defaults to a(1) = 1. Otherwise, filters A and B are of the
% form used in the built-in FILTER function. Filters A and B must have
% the same length (or same number of rows, if they are matrices).
% 
% [Y,Zf] = FILTERINT16(B,A,X,CHANINDEX,Zi) gives access to the initial
% and final conditions. The default for Zi is zeros, corresponding to
% zero-padding the input data. Each _column_ of Zi/Zf specifies the
% conditions for the corresponding row of X. Z is of size
% (filtlen-1)-by-length(CHANINDEX).
%
% Internal calculations are done with float32 precision; the final result is
% rounded off to an int16 to allow natural interface with other int16
% routines.
% Warning: because MATLAB builds filters for float64 precision, this can
% get you into trouble if your filters are too close to the edge of
% instability. Check them out using FREQZ. You can sacrifice speed for
% accuracy by removing the comment on one line of filterint16.c.
%
% See also FILTER, FREQZ.

