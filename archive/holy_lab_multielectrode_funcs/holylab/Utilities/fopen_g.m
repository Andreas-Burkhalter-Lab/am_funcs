function [fid, msg, reused]=fopen_g(varargin)
% fopen_g: a wrapper function to fopen() and test if the fid is reused.
% syntax:
%    [fid, msg, reused]=fopen_g(varargin)
% post:
%    reused: 1 if matlab runs on linux and the file is already open before calling 
%           fopen_g(). (the fid will be same as what last fopen_g() returns)
   oldFids=fopen('all');
   [fid, msg]=fopen(varargin{:});
   reused = ~isempty(find(oldFids==fid));
   