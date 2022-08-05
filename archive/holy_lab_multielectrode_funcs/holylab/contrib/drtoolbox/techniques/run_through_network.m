function result = run_through_network(network, X)
%RUN_THROUGH_NETWORK Runs a set of points through a neural network
%
%   result = run_through_network(network, X)
%
% The function runs a set of points through a neural network that was
% trained by the TRAIN_DEEP_NETWORK function using either the 'Backprop' or
% 'Autoencoder' setting.
%
%
% (C) Laurens van der Maaten
% Maastricht University, 2007
%


    if ~strcmp(network{1}.tied, 'yes')
        disp('This function is meant for networks with undirected weights.');
    end
    
    % Run through network
    no_layers = length(network) + 1;
    activations = X;
    for i=1:no_layers - 1
        if strcmp(network{i}.type, 'sigmoid')
            activations =  1 ./ (1 + exp(-(activations * network{i}.W' + repmat(network{i}.bias', [size(X, 1) 1]))));
        else
            activations = activations * network{i}.W' + repmat(network{i}.bias', [size(X, 1) 1]);
        end
    end
    result = activations;    
    