% Statistics and data analysis
%
% Clustering:
%   adaptive_meanshift - perform locally-adaptive meanshift (LAMS)
%   adaptive_meanshift1 - shift a single landmark once by adaptive meanshift
%   adaptive_meanshift_step - take one step of an adaptive meanshift
%   bn_image_filter  -  balanced neighborhood spatial filtering of images
%   bn_preordered    - the balanced neighborhood criterion
%                      This is the "lowest-level" function, which implements the balanced
%                      neighborhood criterion.
%   bn_preordered_gaussian - balanced neighborhood using Gaussian weights
%   centroid_pls     - find centroid so that power-law scaled data have zero mean
%   choose_landmarks - split a group of points up into local regions
%   climb_gaussian   - cluster by moving points up a mixture-of-gaussians density
%   clust_em_climb   - cluster data based on an EM estimate of the density,
%                      and the climbing of landmarks up to peaks of the density
%   cluster_compare  - compare two clusterings
%   clusterplot      - show clusters as colored dots
%   crossmover       - move a set of representative points to center of mass
%   crossmoverRx     - move landmarks to weighted center of mass
%   find_outliers    - find points too far from landmarks
%   flow_bn_stats    - statistics on convergence locations under balanced-neighborhood mean shift
%   flow_bn_weights  - flow points by balanced neighborhood mean shift, and report weights 
%                      associated with each template
%   iscycling        - test whether a point has entered a cyclic trajectory
%   kcenters         - the k-centers/k-medoids/PAM algorithm
%   kmeans_hard      - fast "hard K-means" algorithm
%   kmeans_relax     - optimize the position of centroids, given a starting guess
%   leapfrog         - rapid convergence for mean shift algorithms by graph iteration
%   make_gaussian_clusters - create simulated data containing clusters
%   maps2sim         - convert a set of MSAMS-maps to a similarity matrix
%   mdcluster        - multidimensional clustering GUI
%   mdcluster_options- set options for mdcluster
%   meanshift        - Use mean-shift algorithm to move landmarks
%   meanshift_errorbars - perform meanshift on a set of points with known uncertainties
%   meanshift1d       - cluster points in one dimension by mean-shift algorithm
%   meanshift1d1      - mex file
%   merge_clust       - manually merge clusters by their cluster #
%   mindist          - Calculate pairs of nearest-neighbors
%   reindex_landmarks - update landmark structure when points are resampled
%   reorder_clust     - re-label clusters, arranging them in a different order
%   reorder_clust_tsp - sort clusters by minimum "travel" distance
%   reorder_manually  - gui that allows you to click on columns and drag them to a 
%                       new position. Click "Done" when satisfied.
%   rmeans            - Fast nonparametric clustering algorithm
%   sim2clust         - convert similarity matrix into a set of clusters
%   vmams_step1       - move a landmark by variable metric adaptive mean shift
%
% Histogramming functions:
%   HistSplit    - Fast histogram (bin boundaries). Also see histc.
%   binpairs     - Binning with potentially-overlapping bins.
%   hist2d       - 2-d histogram.
%   rehist       - Make a coarse histogram from a fine one.
%
% Miscellaneous:
%   ConfInterval - Percentile interval for a parameter.
%   kolsmir      - Kolmogorov-Smirnov statistics for 2 data sets.
%
%   fitstep      - Fit data to a step transition.
%   linregress   - Linear regression.
%   midrank      - Convert values to ranks (for non-parametric stats)
%   connected_components - find the connected components of a graph
%   fit_truncated_gaussian - estimate parameters of Gaussian from incomplete data
%   correct_gaussian_parameters - estimate the complete Gaussian from partial inflation (assumes centered)
%   det_cov      - determinant of the covariance matrix, for different constrained parametrizations
%   dist_mahalanobis - a generalized distance based on a covariance
%   em_gauss     - expectation-maximization (EM) algorithm with spherical gaussians
%   em_gauss_iter - one interation of em algorithm
%   em_relax     -  optimize gaussian mixture model by EM, given a starting guess
%   fdistribution_approx_params - calculate parameters for approximate F distribution
%   find_chisq_boundary - vary fitting parameter until chisq increases by 1
%   fitstep      - fit a sequence of points to a step transition
%   Kolsmir      - kolmogorov-smirnov statistics for 2 data sets
%   LDA          - Linear Discriminant Analysis
%   Linregress   - perform linear regression
%   lognchoosek  - Stirling's approximation for nchoosek
%   metric_biggest_dimensions  - compare points based on largest coordinates
%   medfilt      - Median filtering with correct behavior at data edges.
%   medfilt1     - Median filtering in one dimension.
%   medfilt1improved  - by mathworks  One dimensional median filter
%   moments_from_points - compute zeroth through 2nd moments from a set of points
%   mutual_information - compute the mutual information between labellings
%   PCA          - do Principal Component Analysis and return all components & eigenvalues
%   PCA_absdev   - a variant of PCA robust to outliers
%   PCA_absdev_addone - add a new direction to existing PCA_ABSDEV decomposition
%   PCA_compare  - calculate the degree of overlap between projections
%   pdist2       - Pairwise distance between two sets of observations
%   plsdata      - power-law scale data
%   sqrdist      - find Euclidean square distances between points
%   sqrdist_graph-  square-distances to an off-graph point
%   t2_moments   - compute the moments of T^2 with constrained covariance matrices
%   thres_discrim- compute discriminability as a function of threshold
%   unique_data  - reduce to a set of unique data points & multiplicities
%   varmean      - compute variance and mean & save memory


