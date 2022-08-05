function options = mdcluster_options(input_options)
% MDCLUSTER_OPTIONS: set options for mdcluster
%
% Syntax:
%    options = mdcluster_options
% You can examine the output to learn about the defaults available. Then
% you can modify the fields in "options" to tune to your specific
% application.
%
% The options structure has the following fields:
%   projectionMethods: a structure array of methods (see below)
%   defaultProjectionMethod: the name of the projection method to use
%     upon entry (if missing or blank, no method is used)
%   max_points_for_projection (default Inf): the maximum number of data
%     points to use when calculating new projection directions.
%   clusterMethods: a structure array of methods (see below)
%   blocking (default true): determines whether the GUI blocks (using
%     UIWAIT) while executing.  Setting blocking to true is appropriate
%     if you're calling this from the command line. If it instead is to
%     be embedded in a larger GUI, you can set blocking to be a
%     structure with the following fields:
%       handle:  the handle of a handle graphics object in which to
%         return the final data once completed (the "parent" object)
%       appdata: the name to use for the application data (see
%         SETAPPDATA) in the "parent" object containing the final result
%         (the point labels)
%   showinfo (default true): if true, displays text indicating the # of
%     points, the dimensionality, the current # of clusters.
%   markerSizes (default [6 25]): the sizes of [unselected selected] markers
%
% The structure of a method is:
%   name: the name of the method, used in the pull-down menu text
%   valid: a string (e.g., function name), function handle, or logical
%     used to determine whether the method is currently available. The
%     string 'haveClusters' can be used for a projection method (e.g.,
%     LDA) that needs to have some predefined clusters. Any function
%     specified cannot take arguments, so you may need to use anonymous
%     functions or add functions to the MDEXPLORE internals.
%   execute: the function to call to compute the projected data points
%     (for a projection method) or to calculate the labels of the data
%     points (for a clustering method).  An execute function has the
%     following syntax:
%        output = execute(data,s)
%     where data is a a d-by-N matrix of data points and s is a structure
%     which may hold extra parameters or options for the method.
%   parameterList: a structure array of the form
%     name: contains the text to display in the GUI
%     value: a string representing the value of the parameter
%     castfcn: a function handle used to extract the value from the string
%       Examples of a castfcn:
%         @str2double for an ordinary scalar parameter
%         @(x) sscanf(x,'%d') for a parameter that should be an integer
%           (internally in this function you can use @str2int as a
%           better, and shorter, choice)
%         [] if the parameter is naturally a string (or, if only certain
%           strings are permitted, you can include a function that
%           validates the strings before copying them to the output)
%     The structure s (described above under "execute") is derived from
%     parameterList by executing "castfcn" on "value" and assigning the
%     result to a field "name".
%
% See also: MDCLUSTER.

% Copyright 2008 by Timothy E. Holy

% Note: mdcluster_proj_gateway, drt_gateway, and mdcluster_clust_gateway
% are internal to this function.  

if (nargin < 1)
  input_options = struct;
end
input_options = default(input_options,'use_maxpoints_in_projection',false,'default_n_dims',2,'maxpoints',4000);
default_n_dims = input_options.default_n_dims;

% A simple parameter list entry for the # of dimensions (components)
ndimpL = struct('name','# dimensions','value',default_n_dims,'castfcn',@str2int);
ndimpLmaxpL = ndimpL;
if input_options.use_maxpoints_in_projection
  ndimpLmaxpL(end+1) = struct('name','max points','value',input_options.maxpoints,'castfcn',@str2int);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          projectionMethods            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pM = struct('name','Raw','valid',1,'execute',@(x,s) x,'parameterList',[]);

% For dimensionality reduction toolbox routines, we just append the extra
% parameters. That means that the parameterList has to be created in the
% order expected by dimredtoolbox routines. Start with '# dimensions'.

%PCA
if (exist('pca','file') == 2)
  s = which('pca');
  if ~isempty(strfind(s,'drtoolbox'))
    % Use the DRT method
    pM(end+1) = struct('name','PCA','valid',1,'execute',@(x,s) ...
		       drt_gateway('PCA',x,s),'parameterList',ndimpL);
  else
    % Use our code
    pM(end+1) = struct('name','PCA','valid',1,'execute',@(x,s) ...
		       mdcluster_proj_gateway('PCA',x,s),'parameterList',ndimpLmaxpL);
  end
end

%LDA
if (exist('lda','file') == 2)
  s = which('lda');
  if ~isempty(strfind(s,'drtoolbox'))
    % Use DRT method
    pM(end+1) = struct('name','LDA','valid','haveClusters','execute',@(x,s) ...
		       drt_gateway('LDA',x,s),'parameterList',ndimpL);
  else
    % Use our code (more robust)
    pM(end+1) = struct('name','LDA','valid','haveClusters','execute',@(x,s) ...
		       mdcluster_proj_gateway('LDA',x,s),'parameterList',ndimpLmaxpL);
  end
end

%ICA
if (exist('fastica','file') == 2)
  icapL = ndimpLmaxpL;
  icapL(end+1) = struct('name','nonlinearity','value','tanh','castfcn',@validate_ica_nonlinearity);
  pM(end+1) = struct('name','ICA','valid',1,'execute',@(x,s) ...
		     mdcluster_proj_gateway('ICA',x,s),'parameterList', icapL); %...
% 		     struct('name',{'# dimensions','nonlinearity'}, ...
% 			    'value',{default_n_dims,'tanh'},'castfcn',{@str2int,@validate_ica_nonlinearity}));
end

if (exist('compute_mapping','file') == 2)
  %Isomap
  pM(end+1) = struct('name','Isomap','valid',1,'execute',@(x,s) ...
		     drt_gateway('LandmarkIsomap',x,s),'parameterList', ...
		     struct('name',{'# dimensions','# neighbors','frac. landmarks'}, ...
			    'value',{default_n_dims,12,1},'castfcn',{@str2int,@str2int,@validate_fraction}));
%  pM(end+1) = struct('name','Isomap','valid',1,'execute',@(x,s) ...
%		     drt_gateway('Isomap',x,s),'parameterList', ...
%		     struct('name',{'# dimensions','# neighbors'}, ...
%			    'value',{default_n_dims,12},'castfcn',{@str2int,@str2int}));
  %LPP
  pM(end+1) = struct('name','LPP','valid',1,'execute',@(x,s) ...
		     drt_gateway('LPP',x,s),'parameterList', ...
		     struct('name',{'# dimensions','# neighbors'}, ...
			    'value',{default_n_dims,12},'castfcn',{@str2int,@str2int}));
  %LLC
  pM(end+1) = struct('name','LLC','valid',1,'execute',@(x,s) ...
		     drt_gateway('LLC',x,s),'parameterList', ...
		     struct('name',{'# dimensions','# neighbors','# analyzers'}, ...
			    'value',{default_n_dims,12,40},'castfcn',{@str2int,@str2int,@str2int}));
  %CFA
%  pM(end+1) = struct('name','CFA','valid',1,'execute',@(x,s) ...
%		     drt_gateway('CFA',x,s),'parameterList', ...
%		     struct('name',{'# dimensions','# analyzers'}, ...
%			    'value',{default_n_dims,40},'castfcn',{@str2int,@str2int}));
end



options.projectionMethods = pM;

options.defaultProjectionMethod = 'Raw';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           clusterMethods              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cM = repmat(struct('name','','valid',1,'execute',[],'parameterList',[]),1,0);

% MSAMS
if (exist('msams','file') == 2)
  cM(end+1) = struct('name','MSAMS','valid',1,...
		     'execute', @(x,s) mdcluster_clust_gateway('MSAMS',x,s),...
		     'parameterList',struct('name',{'max points','max landmarks'},'value',{Inf,5000},'castfcn',{@str2int,@str2int}));
end

% k-means
kmpL = struct('name','# clusters','value',2,'castfcn',@str2int);
if (exist('kmeans_hard','file') == 2)
  % Fast custom implementation of k-means
  cM(end+1) = struct('name','k-means','valid',1,...
		     'execute', @(x,s) ...
		     mdcluster_clust_gateway('kmeans_hard',x,s),'parameterList',kmpL);
elseif (exist('kmeans','file') == 2)
  % Use Matlab implementation (slower)
  cM(end+1) = struct('name','k-means','valid',1,...
		     'execute', @(x,s) ...
		     mdcluster_clust_gateway('kmeans_matlab',x,s),'parameterList',kmpL);
end

% GMM
if (exist('gmdistribution','file') == 2)
  cM(end+1) = struct('name','GMM','valid',1,...
		     'execute',@(x,s) ...
		     mdcluster_clust_gateway('GMM',x,s),'parameterList', ...
		     kmpL);
end

% Hierarchical clustering
if (exist('dendrogram','file') == 2)
  cM(end+1) = struct('name','Hierarchical','valid',1,...
		     'execute',@(x,s) ...
		     mdcluster_clust_gateway('hierarch',x,s),'parameterList', ...
		     struct('name','cutoff','value',0.8,'castfcn',@str2double));
end


% Manual merging
cM(end+1) = struct('name','Manual merge','valid',1,...
  'execute',@(x,s) mdcluster_clust_gateway('manualmerge',x,s),...
  'parameterList',[]);

% Manual clustering option
cM(end+1) = struct('name','Manual','valid',1,...
		   'execute',@(x,s) mdclust_manual(x,s),...
		   'parameterList',[]);

options.clusterMethods = cM;


options.blocking = true;
options.showinfo = true;
options.max_points_for_projection = Inf;
options.markerSizes = [6 25];
options.axis_equal = false;

function i = str2int(s)
  i = str2double(s);
  if ~isinf(i)
    i = round(i);
  end
  
function f = validate_fraction(s)
  f = str2double(s);
  if (f > 1)
    f = 1;
  end
  if (f < 0)
    f = 0;
  end
  
function sout = validate_ica_nonlinearity(s)
  sout = lower(s);
  valid_strings = {'pow3','tanh','gauss','skew'};
  if isempty(strmatch(sout,valid_strings,'exact'))
    sout = 'tanh';
  end
  
function output = mdcluster_proj_gateway(methodname,x,s)
  ndimsname = 'n_dimensions';
  ndims = s.(ndimsname);
  ndims = min(ndims,size(x,1));
  if isfield(s,'max_points')
    max_points = s.max_points;
  else
    max_points = Inf;
  end
  switch methodname
   case 'PCA'
    pd = pca(x',struct('maxpoints',max_points));
    pd = pd(:,1:ndims);
    output = pd'*x;
   case 'LDA'
    pd = lda(x,s.labels,struct('maxpergroup',max_points));
    ndims = min(ndims,size(pd,2));
    pd = pd(:,1:ndims);
    output = pd'*x;
   case 'ICA'
    skip = ceil(size(x,2)/max_points);
    if (skip < 1)
      skip = 1;
    end
    if (skip == 1)
      output = fastica(x,'numOfIC',s.(ndimsname),'g',s.nonlinearity);
    else
      [A,W] = fastica(x(:,1:skip:end),'numOfIC',s.(ndimsname),'g',s.nonlinearity);
      output = W*x;
    end
  end
  
  
function output = drt_gateway(methodname,x,s)
% For dimensionality reduction toolbox routines, we just append the extra
% parameters. They therefore have to appear in the parameterList in the
% order described in compute_mapping.
  labels = s.labels;
  s = rmfield(s,'labels');
  c = struct2cell(s)';
  c = [{methodname} c];  % to handle case where s has no fields and c is
                         % therefore empty
  output = compute_mapping(x',c{:})';
  if (size(output,2) < size(x,2))
    warndlg('Not all points were mapped. You might be able to visualize the points, but you may also encounter errors.');
  end

  
function output = mdcluster_clust_gateway(methodname,x,s)
%   % If any are selected, we still want to cluster with awareness of the
%   % whole distribution...
  % Nope, just cluster the selected blob
  if any(s.selected) && ~all(s.selected)
    clabel = agglabel(s.labels);
    selIndex = cat(2,clabel{s.selected});
    x = x(:,selIndex);
  end
  selDone = false;
  switch methodname
   case 'MSAMS'
    skip = ceil(size(x,2)/s.max_points);
    if (skip < 1)
      skip = 1;
    end
    points = x(:,1:skip:end);
    % Choose up to 95% of landmarks randomly, but only when there are many
    % landmarks. Choosing 5% by distance means that we are guaranteed to
    % put some landmarks in the most distant points no matter what.
    lminfo = choose_landmarks(points,s.max_landmarks,...
      struct('frac_random',min(95,s.max_landmarks)/100));
    clust = msams(points,lminfo,struct('flow_only_landmarks',true,'consolidate',false,'reflow',true));
    landmarkClust = clust(lminfo.landmark_xIndex);
    if (skip > 1)
      % For each point, find the closest landmark
      [dist,lmIndex] = mindist(x,lminfo.landmarks);
      % Assign each point to the closest landmark's cluster
      output = landmarkClust(lmIndex);
    else
      output = landmarkClust(lminfo.landmarkAssignment);
    end
   case 'kmeans_hard'
    output = kmeans_hard(x,s.n_clusters,struct('remove_empty',false));
   case 'GMM'
    obj = gmdistribution.fit(x',s.n_clusters,'Regularize',1e-12);
    output = cluster(obj,x');
   case 'hierarch'
    Y = pdist(x');
    Z = linkage(Y);
    figure; dendrogram(Z)
    output = cluster(Z,'cutoff',s.cutoff,'criterion','distance');
   case 'manualmerge'
     selectedIndex = find(s.selected);
     if ~isempty(selectedIndex)
       labels = s.labels;
       labelc = agglabel(s.labels);
       labels(cat(2,labelc{selectedIndex})) = labels(labelc{selectedIndex(1)}(1));
       [ulabels,tmp,output] = unique(labels);
     end
     selDone = true;
  end
  % ... but we apply the results only to the selected cluster(s)
  if any(s.selected) && ~all(s.selected) && ~selDone
    clabel = agglabel(s.labels);
    selIndex = cat(2,clabel{s.selected});
    tmplbl = s.labels;
    %tmplbl(selIndex) = -output(selIndex); % give a guaranteed-different set of labels
    tmplbl(selIndex) = -output; % give a guaranteed-different set of labels
    [ul,tmp,output] = unique(tmplbl);
  end
    
  