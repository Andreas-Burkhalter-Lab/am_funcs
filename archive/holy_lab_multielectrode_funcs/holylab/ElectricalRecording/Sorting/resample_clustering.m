function [clusters mdscaling full_clust_mat] = resample_clustering(indata, options)
% resample_clustering performs clustering using subsets of available observations
% Syntax [clusters cluster_tightness] = resample_clustering(indata, options)
%
% Inputs: 
%         indata (req'd): a n_observationsxn_samples matrix of your data
%         options: struct containing optional subfields
%                .resamp (default min([1 0.5*n_samples])):
%                        =scalar # of observations to be used 1 < resamp < n_observations
%                .iter (default 100): scalar number of iterations larger numbers
%                                     require longer processing time
%                .fxn (default @msams): function to do clustering.  must
%                                       use indata-style input and give
%                                       clust output
%                .display: 'on' or 'off'
% Outputs:
%         clusters: 1xn_samples matrix
%         full_clust_mat: n_samplesxn_samples square comparison matrix
%                         values represent a normalized dissimilarity matrix.
%
% See also msams

% Copyright Julian P. Meeks (Timothy Holy Laboratory)

%% input checking and setup
fixed_sample = false;
input_obs = size(indata,1);
input_samps = size(indata,2);
if ~exist('options', 'var');
    options = struct;
end
options = default(options, 'resamp', max([1 ceil(0.5*input_obs)]));
options = default(options, 'iter', 100);
options = default(options, 'fxn', @msams);
options = default(options, 'fxnops', struct('factor', 1, 'mintocheck', 1));
options = default(options, 'display', 'off');
%options = default(options, 'fxn', @kmeans);
%options = default(options, 'fxnops', struct('factor', 1, 'mintocheck', 1));

if options.resamp < 1 || options.resamp >input_obs
    error('options.resamp is invalid: 1 < options.resamp < size(indata,1)');
else
    resamp = options.resamp;
end

iter = options.iter;
fxn = options.fxn;
fxnops = options.fxnops;
% currently no checks on options.iter or options.fxn
%% 
clust_pairings = zeros(input_samps);
fprintf('Iterating: (of %d)',iter);
if iter > nchoosek(size(indata,1),options.resamp)
    choices = nchoosek(1:size(indata,1),options.resamp);
    iter = size(choices,1);
    fixed_sample = true;
end

for ii = 1:iter
  % select nobs choose resamp manually
  copy = 1:input_obs;
  if fixed_sample == true
      theseobs = choices(ii,:);
  else
      for nki = 1:resamp
          theseobs(nki) = ceil(input_obs*rand);
          while copy(theseobs(nki)) == 0
              theseobs(nki) = ceil(input_obs*rand);
          end
          copy(theseobs(nki))=0;
      end
  end
  % prepare input matrix
  thisdata = indata(theseobs,:);
  % perform clustering on this set
  if strmatch(func2str(fxn), 'msams')
      clust(ii,:) = fxn(thisdata, fxnops);
  elseif strmatch(func2str(fxn), 'kmeans')
      clust(ii,:) = fxn(thisdata', resamp, 'emptyaction', 'drop');
  end
  % add to clust_pairings
  for ci = 1:max(clust(ii,:))
      thisclust = find(clust(ii,:)==ci);
      for tci = 1:size(thisclust,2)
          clust_pairings(thisclust(1:end ~= tci),thisclust(tci)) = clust_pairings(thisclust(1:end ~= tci),thisclust(tci))+1/size(thisclust,2);
      end
  end
  if mod(ii,10) == 0
      fprintf('%d..',ii);
      if mod(ii,100)==0
          fprintf('\n..');
      end
  end
end
fprintf('done!\n');
clust_pairings(find(eye(size(clust_pairings)))) = NaN;
%% use mdscaling to create contrained dimensions upon which to do final clustering
clust_pairings = clust_pairings/max(max(clust_pairings,[],2),[],1); % should now scale ~0-1
clust_pairings = abs(1-clust_pairings); % should invert things linearly
clust_pairings(isnan(clust_pairings))=0;
fprintf('Calculating multidimensional scaling: ');
clust_mdscale = mdscale(clust_pairings,resamp,'options',statset('maxiter',10000), 'start', 'random', 'criterion', 'sstress');
if strmatch(options.display, 'on')
  figure; plot3(clust_mdscale(:,1),clust_mdscale(:,2), clust_mdscale(:,3),'k.');
end
fprintf(' Done!\n');
msams_mdscale = msams(clust_mdscale',struct('factor', 1, 'min_to_check', 3));
%% assign outputs
clusters = msams_mdscale;
mdscaling = clust_mdscale;
full_clust_mat = clust_pairings;
end