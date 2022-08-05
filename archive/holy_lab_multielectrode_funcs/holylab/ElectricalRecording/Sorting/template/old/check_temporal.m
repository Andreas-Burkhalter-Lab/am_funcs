function check_temporal(landmarkClusters, landmarkFinalPos, landmarks, landmarkPos)
% Check the temporal correlation pattern, to see whether the different
% landmarks are just covering different parts of the spike
showtemporal = false;
if showtemporal
  htemporal = figure;
end
hspatial = figure;
for i = 1:length(landmarkClusters)
  ngroups = length(landmarkClusters{i});
  if showtemporal
    figure(htemporal)
    clf
    tc = cell(1,ngroups);
    for j = 1:ngroups
      tc{j} = ts(landmarks{landmarkClusters{i}(j)});
    end
    for j = 1:ngroups
      for k = j:ngroups
        xl = [-5 5];
        npb = crosscorrspike(tc{j},tc{k},xl(2)*10,101);
        subplot(ngroups,ngroups,(j-1)*ngroups+k);
        bar(linspace(xl(1),xl(2),length(npb)),npb);
        set(gca,'XLim',xl);
      end
    end
    %title(num2str(i));
  end % if, showtemporal
  
  figure(hspatial)
  
  % plot cluster_i's landmarks' initial positions
  plot(landmarkPos(:,landmarkClusters{i})); 
  hold on; 
  
  % plot cluster_i's landmarks' final positions (these final positions
  %    should "roughly" same, i.e. converge into one point)
  plot(landmarkFinalPos(:,landmarkClusters{i}),'k'); 
  hold off; 
  
  title(num2str(i));
  numlabels = mat2cell(num2str((1:ngroups)'),ones(ngroups,1),length(num2str(ngroups)));
  legend(numlabels{:})
  pause
end
