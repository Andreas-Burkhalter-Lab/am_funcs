function X = transpose

Z = maxjumpspec2

for i=1:23
    x = Z{i};
    x = x';
    X{i} = x;
    clear x;
end

end
