function network = train_deep_network(X, layers, finetune, targets)
%TRAIN_DEEP_NETWORK Trains a deep (multi-layer) network using RBMs
%
%   network = train_deep_network(X, layers, finetune)
%   network = train_deep_network(X, layers, 'Backprop', targets);
%
% The function trains a deep multi-layer feedforward network on the 
% data specified in X by training each layer separately using Restricted 
% Boltzmann Machine training. The depth and number of nodes in the network 
% are specified by the vector layers. For instance, if layers is set to 
% [100 50 10], a network is trained with three hidden layers with 
% respectively 100, 50 nodes and 10 nodes. The number of input (visual)
% nodes is determined by the dimensionality of the input data X. The
% network is trained using a greedy approach. The network may be finetuned
% using backpropagation or contrastive wake-sleep by setting finetune. 
% Possible values are 'Backprop', 'WakeSleep', or 'None' (default = 'None').
% The network is returned in the cell-array network.
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.4b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007


    if ~exist('layers', 'var') || isempty(layers)
        layers = [20 10];
    end
    if ~exist('finetune', 'var') || isempty(finetune)
        finetune = 'None';
    end
    disp(' ');
    
    % Initialize some variables
    no_layers = length(layers);
    network = cell(1, no_layers);
    X = X -  min(min(X));
    X = X ./ max(max(X));
    origX = X;
            
    % Learn layer-by-layer to get an initial network configuration
    for i=1:no_layers
        
        disp(['Training layer ' num2str(i) '...']);
        
        % Train current layer
        if i ~= no_layers
            network{i} = train_rbm(X, layers(i), 'sigmoid');
        else
            if ~strcmp(finetune, 'Backprop')
                network{i} = train_rbm(X, layers(i), 'linear');
            else
                network{i}.W = randn([layers(i - 1) layers(i)]) * 0.1;
                network{i}.bias_upW = zeros([1 layers(i)]);
                network{i}.bias_downW = zeros([1 layers(i - 1)]);
                network{i}.type = 'linear';
            end
        end
        
        % Transform data using learned weights
        if i ~= no_layers
            X = 1 ./ (1 + exp(-(X * network{i}.W + repmat(network{i}.bias_upW, [size(X, 1) 1]))));
        end
    end

    % Perform finetuning if desired
    switch finetune
        
        case 'WakeSleep'
            disp('Finetuning the network using the contrastive wake-sleep algorithm...');
            network = wake_sleep(network, origX); 
            
        case 'Autoencoder'
            disp('Finetuning the autoencoder using backpropagation...');
            no_layers = length(network);
            for i=1:no_layers
                network{2 * no_layers + 1 - i} = network{i};
                network{i}.W = network{i}.W';
                network{i}.bias = network{i}.bias_upW';
                network{2 * no_layers + 1 - i}.bias = network{i}.bias_downW';
            end
            network{no_layers + 1}.type = 'sigmoid';
            network = backprop(network, origX, origX);
            
        case 'Backprop'
            disp('Finetuning the network using backpropagation...');
            for i=1:length(network)
                network{i}.W = network{i}.W';
                network{i}.bias = network{i}.bias_upW';
            end
            network = backprop(network, origX, targets);
            
        case 'Unsupervised'
            for i=1:length(network)
                network{i}.W = network{i}.W';
                network{i}.bias = network{i}.bias_upW';
            end
            
        otherwise
            % Do nothing
    end
