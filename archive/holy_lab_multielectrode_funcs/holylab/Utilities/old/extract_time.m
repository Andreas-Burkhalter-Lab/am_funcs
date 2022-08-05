function [totaltime, simplifiedTimeStr]=extract_time(a_time_string)
% extract_time: parse a string for time
%                   
% Syntax:
%    [totaltime, simplifiedTimeStr]=extract_time(a_time_string)
%    
% @pre:   
%    a_time_string: a one-row string.
%  
% @post:
%    totaltime: total time in seconds
%    simplifiedTimeStr: a string only has time info
%    
% @eg:
%    [totaltime, simplifiedTimeStr]=extract_time('30s 10m and 2 min');
%    we got: totaltime=750, simplifiedTimeStr='30s 10m 2 min '
% 
% @see: 
%    stand-alone executable: extract_time

