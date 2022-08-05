function Y = sym_sne(X, d)
%SNE Implementation of symmetric Stochastic Neighbor Embedding
%
%   Y = sym_sne(X, d)
%
% Runs the symmetric Stochastic Neighbor Embedding algorithm. The 
% high-dimensional datapoints are specified by X. The target dimensionality
% if specified in d. The function returns the embedded points in Y.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.4b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

    % Initialize some variables
    N = size(X, 1);                 % number of instances
    sigma = 1;                      % variance of Gaussian kernel
    eta = .5;                       % learning rate
    max_iter = 4000;                % maximum number of iterations

    % Normalize dataset
    X = X - min(min(X));
    X = X / max(max(X));
    
    % Initialize embedding coordinates randomly (close to origin)
    Y = rand(N, d) * .01;
    dJ = zeros(size(Y, 1), d);
    
    % Probability computations for the high-dimensional space
    P = exp(-L2_distance(X', X') .^ 2) ./ (2 * sigma ^ 2));
    P = P / (sum(sum(tril(P))) - sum(diag(P)));

    % Iterating loop
    for iter=1:max_iter
        
        % Compute gradient of J
        for i=1:N
            dJ(i,:) = sum((repmat(Y(i,:), N, 1) - Y) .* repmat(P(i,:)', 1, d), 1);        
        end
     
        % Steepest descent
        Y = Y - eta * dJ;
		
		% Normalize Y (in order to prevent very large coordinates)
		Y = Y / max(max(Y));
        
        % Add some Gaussian noise (reduce jitter over time)
        Y = Y + (randn(size(Y)) * (.3 - iter * (.3 / max_iter)));
        
        % Display information
        disp(['Iteration ' num2str(iter) '...']);
    end
    