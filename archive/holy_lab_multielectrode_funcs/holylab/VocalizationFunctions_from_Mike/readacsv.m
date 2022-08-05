function csvcell = readacsv(filename)

x = textread(filename,'%s','delimiter','\n');
x = x';
row = cell(1,length(x));
lengths = zeros(1,length(x));

    for i = 1:length(x)
        row{i} = regexp(x{i},',','split');
        lengths(i) = length(row{i});
    end
    
maxlength = max(lengths);
csvcell = cell(length(row),maxlength);

for i = 1:size(csvcell,1)
    for j = 1:lengths(i)
    csvcell{i,j} = row{i}{j};
    end
end



end
