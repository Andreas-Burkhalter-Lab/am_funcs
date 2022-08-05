p(1) = struct('r1',20,'r2',1,'t1',0.020,'t2',0.050);
%p(2) = struct('r1',80,'r2',5,'t1',0.015,'t2',0.035);
p(2) = struct('r1',15,'r2',1.5,'t1',0.015,'t2',0.035);
p(3) = struct('r1',10,'r2',0.5,'t1',0.020,'t2',0.060);
p(4) = p(1);
p(5) = p(3);
minspikes = 5000;
lp = length(p);
for i = 1:lp
  mt(i) = spikesim(p(i));
end
maxtime = max(mt);
n = round(maxtime*minspikes./mt);
for i = 1:lp
  tt{i} = spikesim(p(i),n(i));
  t{i} = cumsum(tt{i});
  tend(i) = t{i}(end);
end
% Insure that the spike trains cover the same period of time
mintend = min(tend);
for i = 1:lp
  t{i} = t{i}(find(t{i} < mintend));
end
% Now that we've got the spike trains of individual cells, let's make some dirty clusters
% Cluster 1: 80% of the spikes of cell 1.
% Cluster 2: 20% of the spikes of cell 1, and all of the spikes from cells 2 and 3.
% Cluster 3: 80% of the spikes from cell 4, 20% of the spikes from cell 5
% Cluster 4: 20% of the spikes from cell 4, 80% of the spikes from cell 5
contmtrx = [.8 0 0 0 0; .2 1 1 0 0; 0 0 0 .8 .2; 0 0 0 .2 .8];
indx = randperm(length(t{1}));
splitpt = round(0.8 * length(indx));
indx1 = sort(indx(1:splitpt));
indx2 = sort(indx(splitpt+1:end));
c1 = t{1}(indx1);
c2 = sort([t{1}(indx2),t{2},t{3}]);
indx = randperm(length(t{4}));
splitpt = round(0.8 * length(indx));
indx1 = sort(indx(1:splitpt));
indx2 = sort(indx(splitpt+1:end));
c3 = t{4}(indx1);
c4 = t{4}(indx2);
indx = randperm(length(t{5}));
splitpt = round(0.2 * length(indx));
indx1 = sort(indx(1:splitpt));
indx2 = sort(indx(splitpt+1:end));
c3 = sort([c3,t{5}(indx1)]);
c4 = sort([c4,t{5}(indx2)]);
allcs = {c1,c2,c3,c4};

totfrac = zeros(4,5);
for i = 1:4
  for j = 1:5
    totfrac(i,j) = length(t{j})/length(allcs{i});
  end
end
totfrac = totfrac.*contmtrx;
totfrac

% Now calculate the omega parameters
mixparams = struct('tref',.015,'tp',[0.15 0.3],'plot',1);
omega = spikemixing(allcs,mixparams)