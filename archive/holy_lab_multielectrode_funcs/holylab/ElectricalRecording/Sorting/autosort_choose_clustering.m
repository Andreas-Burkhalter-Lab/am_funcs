% This doesn't do anything useful, just for testing purposes.

% Apply landmark cluster numbers to the individual points (so we can
% compare across replicates)
plabel = cell(1,sr.nRfactors);
for i = 1:sr.nRfactors
  for j = 1:sr.nreplicates
    plabel{i}(j,:) = sr.landmarkClust{j}(i,sr.landmarkIndex(j,:));
  end
end
% Calculate the mutual information between all pairs of replicates with the
% same Rfactor
mi = cell(1,sr.nRfactors);
for k = 1:sr.nRfactors
  mi{k} = zeros(sr.nreplicates);
  for i = 1:sr.nreplicates
    for j = i:sr.nreplicates
      mi{k}(i,j) = mutual_information(plabel{k}([i j],:));
      mi{k}(j,i) = mi{k}(i,j);
    end
  end
end
% Calculate the mean mutual information as a function of Rfactor, excluding
% the self mutual information
for k = 1:sr.nRfactors
  tmp_selfmi = diag(mi{k});
  tmp_mi = mi{k} - diag(tmp_selfmi);
  meanmi(k) = mean(tmp_mi(:));
  selfmi(k) = mean(tmp_selfmi);
end

% Calculate the mutual information within a replicate but across labellings
% with different Rfactors
miRf = zeros([sr.nRfactors sr.nRfactors sr.nreplicates]);
for k = 1:sr.nreplicates
  for i = 1:sr.nRfactors
    for j = i:sr.nRfactors
      miRf(i,j,k) = mutual_information([plabel{i}(k,:); plabel{j}(k,:)]);
      miRf(j,i,k) = miRf(i,j,k);
    end
  end
end
% Calculate the mean mutual information as a function of Rfactor
meanmiRf = mean(miRf,3);
