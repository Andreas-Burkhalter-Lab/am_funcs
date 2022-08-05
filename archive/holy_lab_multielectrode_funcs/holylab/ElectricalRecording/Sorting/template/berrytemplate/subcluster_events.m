function [subclusters, aborted]=subcluster_events(snip, separationPositions, options)
% PRE:
%    snip: each column is a snippet; each snippet is saved channel-by-channel
% POST: 
%    subclusters: a cell array of column indices to snip
   
   % pca 
   [pd,sv] = pca(snip'); % NOTE: pd: project directions, each col is a principal components
   proj = pd(:,1:2)'*snip; % NOTE: project snippets on the two most significant principal components
			   % NOTE: proj(1,:) is the project of all snippets on the 1st priniciple component
   %figure; plot(proj(1,:),proj(2,:),'x')
   
   figPC=figure; plot(pd(:,1:2)) % also plot the principal components 
   set(figPC, 'name', 'principal components');
   if isstruct(options.win_locations)
      set(figPC,'Position',options.win_locations.PC)
   end

   % interactively check snippets
   for i = 1:size(snip,2); 
      snipcell{i} = snip(:,i); 
   end
   subClustOps = struct('separationPositions', separationPositions);
   subClustOps.win_locations = options.win_locations;
   figMdexplore = figure;
   subClustOps.fignum = figMdexplore;
   if isstruct(options.win_locations)
       set(figMdexplore,'Position',options.win_locations.mdexplore)
   end
   [subclusters, aborted]=mdexplore2(proj(1:2,:),snipcell, subClustOps);
   %figMdexplore=gcf;
   
   if ishandle(figMdexplore)
     close(figMdexplore); % will trigger the closeRequestFcn
   end
   
   free([figPC]);
   