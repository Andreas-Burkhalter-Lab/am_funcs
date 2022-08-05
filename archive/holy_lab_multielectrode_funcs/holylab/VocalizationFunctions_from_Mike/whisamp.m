function [y,k] = whisamp(snips,nfreq)

n = length(snips);
y = zeros(1,nfreq);

    for i = 1:n
    
    mag = sum(full(abs(snips{i})),2);
    y = y+mag';
    
    end

k = find(y);
y = y(k);

end
