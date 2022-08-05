function network = backprop(network, X, T, max_iter)
%BACKPROP Trains a network on a dataset using backpropagation
%
%   network = backprop(network, X, T, max_iter)
%
% The function trains the specified network using backpropagation on
% dataset X with targets T for max_iter iterations. The dataset X is an NxD
% matrix, whereas the targets matrix T has size NxM. The variable network
% is a cell array that may be obtained from the TRAIN_DEEP_NETWORK
% function. The function returns the trained network in network.

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.4b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007

    if ~exist('max_iter', 'var') || isempty(max_iter)
        max_iter = 100;
    end

    % Initialize some variables
    n = size(X, 1);
    no_layers = length(network);
    batch_size = 1 + round(n / 20);
    
    % Perform the backpropagation
    for iter=1:max_iter
        
        % Print progress
        if rem(iter, 10) == 0
            disp(['Iteration ' num2str(iter) '...']);
        end
        
        % Loop over all batches
        index = randperm(n);
        for batch=1:batch_size:n
            
            % Select current batch
            tmpX = X(index(batch:min([batch + batch_size - 1 n])),:);
            tmpT = T(index(batch:min([batch + batch_size - 1 n])),:);
            
            % Convert the weights and store them in the network
            v = [];
            for i=1:length(network)
                v = [v; network{i}.W(:); network{i}.bias(:)];
            end
            
            % Conjugate gradient minimization using 3 linesearches
            v = minimize(v, 'backprop_gradient', 3, network, tmpX, tmpT);
            
            % Deconvert the weights and store them in the network
            ind = 1;
            for i=1:no_layers
                network{i}.W    = reshape(v(ind:ind - 1 + numel(network{i}.W)),    size(network{i}.W));     ind = ind + numel(network{i}.W);
                network{i}.bias = reshape(v(ind:ind - 1 + numel(network{i}.bias)), size(network{i}.bias));  ind = ind + numel(network{i}.bias);
            end
        end
    end
    