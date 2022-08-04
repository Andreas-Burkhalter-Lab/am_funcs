%% initialize the clustering with pixel assigments

Upix = reshape(Upix, Ly, Lx, nPC);

% we'll come to this later (high-pass filtering to reduce effect of
% neuropil)
for j = 1:nPC
   Upix(:,:,j)  = Upix(:,:,j) - my_conv2(Upix(:,:,j), 10, [1 2]);
end
%% clustering with a fixed number of clusters (k-means + scaling variables)
% decide how many clusters
Nclusters = 1000;

figure('Position', [100 100 800 700])

[Ly, Lx, NT] = size(mov);

% use a Suite2p function for initializing clusters with squares
iclust = initialize_clusters(Upix, Nclusters, 'squares', Lx, Ly); 

Upix = reshape(Upix, [], nPC);

% r is a random hue assigned to each cluster for visualization
r = rand(Nclusters,1);
% saturation is always one for visualization
Sat = ones(Ly, Lx);

% run niter iterations of the algorithm 
niter = 10;
for iter = 1:niter
    lam = zeros(nPC, Nclusters);
     for j = 1:Nclusters
        % these are all the pixels in a cluster
        ipix = (iclust==j);
        
        % determine the mean activity of this cluster
        lam(:, j) = mean(Upix(ipix, :), 1);
     end
    
    % normalize the cluster means to unit norm
    lam = normc(lam);
    
    % compute the product between each cluster mean and each pixel
    cc = Upix * lam;
    
    % reassign pixels to clusters they are most correlated to
    [cmax, iclust] = max(cc, [], 2);
    
    % show pseudo-color figure
    V = max(0, min(.5 * reshape(cmax, Ly, Lx)/mean(cmax(:)), 1));
    H = reshape(r(iclust), Ly, Lx);
    rgb_image = hsv2rgb(cat(3, H, Sat, V));
    imagesc(rgb_image)
    axis off
    drawnow
end

%% how to determine the number of cells?
% let's compute a so-called "correlation map"
% residual is smoothed at every iteration
Upix = reshape(Upix, Ly, Lx, nPC);

% smoothing scale (in pixels)
sig = 2;

us = my_conv2_circ(Upix, sig, [1 2]);

% compute total activity power at each location
V = sq(mean(us.^2,3));
V = double(V);

um = sq(mean(Upix.^2,3));
um = my_conv2_circ(um, sig, [1 2]);

V = V./um ;
%     V = log(V./um);
V = double(V);

figure; imagesc(V); title('correlation map');
%%

