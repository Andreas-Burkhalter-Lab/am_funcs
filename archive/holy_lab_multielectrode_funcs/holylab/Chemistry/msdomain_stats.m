function domain_stats = msdomain_stats(Ms,domains)
  n_domains = length(domains);
  tmp = struct('intensity',0,'mean',[0;0],'std',[0;0],'range',[0;0]);
  domain_stats(n_domains) = tmp;
  c = cell(1,2);
  for i = 1:n_domains
    [c{:}] = ind2sub(size(Ms),domains{i});
    thisMs = Ms(domains{i});
    tmp.intensity = sum(thisMs);
    tmp.std = zeros(2,1);
    for j = 1:2
      tmp.mean(j) = sum(thisMs .* c{j})/tmp.intensity;
      dc = c{j} - tmp.mean(j);
      tmp.std(j) = sqrt(sum(thisMs .* dc.^2)/tmp.intensity);
    end
    tmp.range(1) = range(c{1});
    tmp.range(2) = range(c{2});
    domain_stats(i) = tmp;
  end
end
