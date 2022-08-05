function [s, e] = strip_se(se_in)
% STRIP_SE parses s%e% notation to return the %values
% 
% Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)
%

if iscell(se_in)
    n_se = size(se_in,2);
    s = zeros(size(se_in));
    e = ones(size(se_in));
else
    se_in = {se_in};
    n_se = 1;
    s = 0;
    e = 1;
end

for idx = 1:n_se
    sidx = strfind(se_in{idx},'s');
    counter = 1;
    while ~isempty(str2num(se_in{idx}(sidx+counter)))
        sstr(counter) = se_in{idx}(sidx+counter);
        counter = counter+1;
    end
    s(idx) = str2num(sstr);
    sstr = char;
    eidx = strfind(se_in{idx},'e');
    counter = 1;
    while ~isempty(str2num(se_in{idx}(eidx+counter)))
        estr(counter) = se_in{idx}(eidx+counter);
        if eidx+counter+1 <= size(se_in{idx},2)
            counter = counter+1;
        else
            break
        end
    end
    e(idx) = str2num(estr);
    estr = char;
end

% sort by s, then by e
sortmat = sortrows([s' e']);
col1 = sortmat(:,1);
col2 = sortmat(:,2);
s = col1';
e = col2';
    

end