%%%%%%%%%%%%%%%5
%% This script builds the network described in Kinouchi and Copelli (2006) and sets initial parameters. It should be called by other programs. 
%%%%%%%%%%%%%%
clc; 
clear;
close all;

N = 500;                    % number of elements in the network
sigma=1                     % average local branching ratio
k=10;                        % average connectivity
pmax=2*sigma/k;
states =7;                   % number of refractory states minus 3
T = 1000;                    % duration of the trial for trials with background noise
T1 = 1000;                  % time to run background noise before applying external stimulus dr
r = 0;                      % general value for external input - can take on values of r0 (background noise) and be increased by dr
r0 = 0;                         % original background noise.... average number of firings due to external input per timestep
tau = 750;                     % time constant for recovery of dynamic synapses
u0mean = 0.07;                     % short-term depression parameter from Levina et al. 2007 - fraction of current weight subtracted at spike

v=zeros(1,(N-1)*N/2);                % this vector will be used to form random connections between elements
v(1,randperm((N-1)*N/2,k*N/2))=1;
v = pmax * v .* (rand(1,length(v))); % each connection has random strength from 0 to pmax

W=zeros(N);                          % create weight matrix
W((triu(ones(N))-eye(N))==1)=v;  
W=W+W';                              % copy all connections over the diagonal to make them symmetrical
W = W.*sigma/mean(sum(W));          % multiplies all connections in the network by a factor to ensure that mean(sum(W)) starts as sigma 
W0 = W;

pause on