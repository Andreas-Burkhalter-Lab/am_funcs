function svdtemplates(originalFineClustFilename)
% SVDTEMPLATES: find an indepedent set of template components
% Syntax:
%   svdtemplates(fineClusterFilename)
% 
% Here, fineClusterFilename is a string containing the name of a
% *.fine_cluster file, i.e., the output of cluster_events_fine.
%
% This will save a file with the same basename with extension
% file_cluster_svd.
%
% See also: FIT_COMPONENTS, CLUSTER_EVENTS_FINE.
  
% Copyright 2007 by Timothy E. Holy
  
  load('-mat',originalFineClustFilename);
  
  tpl = cat(2,templates{:});
  % Create a smoothed version of the templates, so that noise doesn't
  % dominate (noise comes mostly from with very few spike waveform
  % examples)
  svd_smoothed = false;
  if svd_smoothed
    [b,a] = butter(2,0.2,'low');
    tpl_sm = filtfilt(b,a,tpl);
    % Do SVD on the smoothed templates
    [U,S,V] = svd(tpl_sm,'econ');
  else
    [U,S,V] = svd(tpl,'econ');
  end
  % Show the user the singular values
  figure
  s = diag(S);
  plot(s,'x');
  % Get input from the user about the number of components to use
  n_components = input('How many components should I keep? ');
  % Find the projections of the kept components
  projections = S(1:n_components,1:n_components)*V(:,1:n_components)';
  if svd_smoothed
    % Now find the directions in terms of the _original_ (unsmoothed)
    % templates that best recreate this set of projections
    W = (tpl*projections')/(projections*projections');
    % Finally, calculate the projections of original templates on these
    % directions (should be very close to first projections)
    projections = W'*tpl;
  else
    W = U(:,1:n_components);
  end
    
  % Define these new directions as the templates, in the format expected
  % for a fine_cluster file.
  templates = cell(1,n_components);
  for i = 1:n_components
    templates{i} = W(:,i);
  end
  
  % Save the results
  [pathstr,basename] = fileparts(originalFineClustFilename);
  newName = fullfile(pathstr,[basename '.fine_cluster_svd']);
  save('-mat',newName,...
       'channels',...
       'fineClusterOptions',...
       'originalFineClustFilename',...
       'medv',...
       'rawClusterFilename',...
       'templates',...
       'thresh',...
       'projections',...
       'sniprange');
  
  