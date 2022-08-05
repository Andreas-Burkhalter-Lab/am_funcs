function [C, dC] = backprop_gradient(v, network, X, targets)
%BACKPROP Compute the cost gradient for CG optimization of a neural network
%
%   [C, dC] = backprop_gradient(v, network, X, targets)
%
% Compute the value of the cost function, as well as the corresponding 
% gradient for conjugate-gradient optimization of a neural network.

% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.4b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007


    % Initialize some variables
    n = size(X, 1);
    no_layers = length(network);

    % Deconvert the weights and store them in the network
    ind = 1;
    for i=1:no_layers
        network{i}.W    = reshape(v(ind:ind - 1 + numel(network{i}.W)),    size(network{i}.W));     ind = ind + numel(network{i}.W);
        network{i}.bias = reshape(v(ind:ind - 1 + numel(network{i}.bias)), size(network{i}.bias));  ind = ind + numel(network{i}.bias);
    end
    
    % Run the data through the network
    activations = cell(1, no_layers + 1);
    activations{1} = [X ones(n, 1)];
    for i=1:no_layers
        if strcmp(network{i}.type, 'sigmoid')
            activations{i + 1} = 1 ./ (1 + exp(-(activations{i} * [network{i}.W'; network{i}.bias'])));
        else
            activations{i + 1} = activations{i} * [network{i}.W'; network{i}.bias'];
        end
        activations{i + 1} = [activations{i + 1} ones(n, 1)];
    end 

    % Compute value of cost function (= cross entropy)
    C = (-1 / n) .* sum(sum(targets  .* log( activations{end}(:,1:end - 1)) + ...
                       (1 - targets) .* log(-activations{end}(:,1:end - 1) + 1)));
    
    % Compute gradients 
    dW = cell(1, no_layers);
    db = cell(1, no_layers);
    Ix = (activations{end}(:,1:end - 1) - targets) ./ n;            % discrepancy between activations and targets (assumes linear output layer)
    for i=no_layers:-1:1
        
        % Compute erros in top layer
        delta = activations{i}' * Ix;
        dW{i} = delta(1:end - 1,:)';
        db{i} = delta(end,:)';
        
        % Backpropagate error
        if i > 1
            if strcmp(network{i - 1}.type, 'sigmoid') 
                Ix = (Ix * [network{i}.W network{i}.bias]) .* activations{i} .* (1 - activations{i});
            else
                Ix = Ix * [network{i}.W network{i}.bias];
            end
            Ix = Ix(:,1:end - 1);
        end
    end

    % Convert gradient information
    dC = repmat(0, [numel(v) 1]);
    ind = 1;
    for i=1:no_layers
        dC(ind:ind - 1 + numel(dW{i})) = dW{i}(:); ind = ind + numel(dW{i});
        dC(ind:ind - 1 + numel(db{i})) = db{i}(:); ind = ind + numel(db{i});
    end
