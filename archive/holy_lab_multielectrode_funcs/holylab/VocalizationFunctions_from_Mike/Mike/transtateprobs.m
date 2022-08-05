function [observed,expected,transitions] = transtateprobs(x)

% Count number of categories
cats = min(x):max(x);
observed = zeros(numel(cats)*numel(cats),1);
expected = zeros(numel(cats)*numel(cats),1);
transitions = zeros(numel(cats)*numel(cats),2);
c = 1;


for i = 1:numel(cats)
    for j = 1:numel(cats)
        observed(c) = numel(find(x(1:end-1) == cats(i) & x(2:end) == cats(j)));
        expected(c) = (numel(find(x == cats(i)))/length(x))*numel((find(x == cats(j)))/length(x))*length(x);
        transitions(c,1) = cats(i);
        transitions(c,2) = cats(j);
        c = c+1;
    end
end

end