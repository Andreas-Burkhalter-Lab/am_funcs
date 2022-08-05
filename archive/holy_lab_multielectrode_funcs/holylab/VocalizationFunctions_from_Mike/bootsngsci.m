function [mns,cis,bootdist] = bootsngsci(x,N)


n = numel(x);


for j = 1:length(x{1})
    for i=1:n
      if ~isempty(x{i})
      y(i)=x{i}(j);
      end
    end
    y = y(~isinf(y));
    [ci,bootdist{j}] = bootci(N,@mean,y);
    mn = mean(y(~isinf(y)));
    
    cis{j}=ci;
    mns(j)=mn;
    clear mn;
    clear ci;
    clear y;
end


end
