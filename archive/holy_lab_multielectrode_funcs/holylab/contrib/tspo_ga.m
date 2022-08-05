function varargout = tspo_ga(xy,dmat,pop_size,num_iter,show_prog,show_res,penalty_func)
%TSPO_GA Open Traveling Salesman Problem (TSP) Genetic Algorithm (GA)
%   Finds a (near) optimal solution to a variation of the TSP by setting up
%   a GA to search for the shortest route (least distance for the salesman
%   to travel to each city exactly once without returning to the starting city)
%
% Summary:
%     1. A single salesman travels to each of the cities but does not close
%        the loop by returning to the city he started from
%     2. Each city is visited by the salesman exactly once
%
% Input:
%     XY (float) is an Nx2 matrix of city locations, where N is the number of cities
%     DMAT (float) is an NxN matrix of point to point distances/costs
%     POP_SIZE (scalar integer) is the size of the population (should be divisible by 4)
%     NUM_ITER (scalar integer) is the number of desired iterations for the algorithm to run
%     SHOW_PROG (scalar logical) shows the GA progress if true
%     SHOW_RES (scalar logical) shows the GA results if true
%     PENALTY_FUNC is an optional penalty function of the syntax
%           d = penalty_func(sequence)
%       that gets added to the distance penalty; sequence is the route, an
%       ordering of 1:n_cities.
%
% Output:
%     OPT_RTE (integer array) is the best route found by the algorithm
%     MIN_DIST (scalar float) is the cost of the best route
%
% 2D Example:
%     n = 50;
%     xy = 10*rand(n,2);
%     pop_size = 60;
%     num_iter = 1e4;
%     show_prog = 1;
%     show_res = 1;
%     a = meshgrid(1:n);
%     dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),n,n);
%     [opt_rte,min_dist] = tspo_ga(xy,dmat,pop_size,num_iter,show_prog,show_res);
%
% 3D Example:
%     n = 50;
%     xyz = 10*rand(n,3);
%     pop_size = 60;
%     num_iter = 1e4;
%     show_prog = 1;
%     show_res = 1;
%     a = meshgrid(1:n);
%     dmat = reshape(sqrt(sum((xyz(a,:)-xyz(a',:)).^2,2)),n,n);
%     [opt_rte,min_dist] = tspo_ga(xyz,dmat,pop_size,num_iter,show_prog,show_res);
%
% See also: tsp_ga, tsp_nn, tspof_ga, tspofs_ga, distmat
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.2
% Release Date: 6/2/09
%
% Penalty function added by Timothy E. Holy

% Process Inputs and Initialize Defaults
nargs = 6;
for k = nargin:nargs-1
    switch k
        case 0
            xy = 10*rand(50,2);
        case 1
            N = size(xy,1);
            a = meshgrid(1:N);
            dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),N,N);
        case 2
            pop_size = 100;
        case 3
            num_iter = 1e4;
        case 4
            show_prog = 1;
        case 5
            show_res = 1;
        otherwise
    end
end

% Verify Inputs
[N,dims] = size(xy);
[nr,nc] = size(dmat);
if N ~= nr || N ~= nc
    error('Invalid XY or DMAT inputs!')
end
n = N;
if (nargin < 7)
  penalty_func = @(order) 0;
end

% Sanity Checks
pop_size = 4*ceil(pop_size/4);
num_iter = max(1,round(real(num_iter(1))));
show_prog = logical(show_prog(1));
show_res = logical(show_res(1));

% Initialize the Population
pop = zeros(pop_size,n);
for k = 1:pop_size
    pop(k,:) = randperm(n);
end

% Run the GA
global_min = Inf;
total_dist = zeros(1,pop_size);
dist_history = zeros(1,num_iter);
tmp_pop = zeros(4,n);
new_pop = zeros(pop_size,n);
if show_prog
    pfig = figure('Name','TSPO_GA | Current Best Solution','Numbertitle','off');
end
for iter = 1:num_iter
    % Evaluate Each Population Member (Calculate Total Distance)
    for p = 1:pop_size
        d = 0; % Open Path
        for k = 2:n
            d = d + dmat(pop(p,k-1),pop(p,k));
        end
        total_dist(p) = d + penalty_func(pop(p,:));
    end

    % Find the Best Route in the Population
    [min_dist,index] = min(total_dist);
    dist_history(iter) = min_dist;
    if min_dist < global_min
        global_min = min_dist;
        opt_rte = pop(index,:);
        if show_prog
            % Plot the Best Route
            figure(pfig);
            if dims == 3, plot3(xy(opt_rte,1),xy(opt_rte,2),xy(opt_rte,3),'r.-');
            else plot(xy(opt_rte,1),xy(opt_rte,2),'r.-'); end
            title(sprintf('Total Distance = %1.4f, Iteration = %d',min_dist,iter));
        end
    end

    % Genetic Algorithm Operators
    rand_pair = randperm(pop_size);
    for p = 4:4:pop_size
        rtes = pop(rand_pair(p-3:p),:);
        dists = total_dist(rand_pair(p-3:p));
        [ignore,idx] = min(dists);
        best_of_4_rte = rtes(idx,:);
        ins_pts = sort(ceil(n*rand(1,2)));
        I = ins_pts(1);
        J = ins_pts(2);
        for k = 1:4 % Mutate the Best to get Three New Routes
            tmp_pop(k,:) = best_of_4_rte;
            switch k
                case 2 % Flip
                    tmp_pop(k,I:J) = fliplr(tmp_pop(k,I:J));
                case 3 % Swap
                    tmp_pop(k,[I J]) = tmp_pop(k,[J I]);
                case 4 % Slide
                    tmp_pop(k,I:J) = tmp_pop(k,[I+1:J I]);
                otherwise % Do Nothing
            end
        end
        new_pop(p-3:p,:) = tmp_pop;
    end
    pop = new_pop;
end

if show_res
    % Plots the GA Results
    figure('Name','TSPO_GA | Results','Numbertitle','off');
    subplot(2,2,1);
    if dims == 3, plot3(xy(:,1),xy(:,2),xy(:,3),'k.');
    else plot(xy(:,1),xy(:,2),'k.'); end
    title('City Locations');
    subplot(2,2,2);
    imagesc(dmat(opt_rte,opt_rte));
    title('Distance Matrix');
    subplot(2,2,3);
    if dims == 3, plot3(xy(opt_rte,1),xy(opt_rte,2),xy(opt_rte,3),'r.-');
    else plot(xy(opt_rte,1),xy(opt_rte,2),'r.-'); end
    title(sprintf('Total Distance = %1.4f',min_dist));
    subplot(2,2,4);
    plot(dist_history,'b','LineWidth',2);
    title('Best Solution History');
    set(gca,'XLim',[0 num_iter+1],'YLim',[0 1.1*max([1 dist_history])]);
end

% Return Outputs
if nargout
    varargout{1} = opt_rte;
    varargout{2} = min_dist;
end
