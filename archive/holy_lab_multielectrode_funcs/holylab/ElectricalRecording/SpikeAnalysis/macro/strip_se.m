function [s, e] = strip_se(se_in)
% STRIP_SE parses s%e% notation to return the %values
% 
% Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)
%

if iscell(se_in)
    n_se = size(se_in,2);
    s = zeros(size(se_in));
    e = zeros(size(se_in));
else
    se_in = {se_in};
    n_se = 1;
    s = 0;
    e = 0;
end

for idx = 1:n_se
    sidx = strfind(se_in{idx},'s');
    counter = 1;
    while ~isempty(str2num(se_in{idx}(sidx+counter))) || ~isempty(strmatch(se_in{idx}(sidx+counter),'.', 'exact')) || ~isempty(strmatch(se_in{idx}(sidx+counter),'_', 'exact'))
        if se_in{idx}(sidx+counter) == '_'
            sstr(counter) = '.';
        else
            sstr(counter) = se_in{idx}(sidx+counter);
        end
        counter = counter+1;
    end
    s(idx) = str2num(sstr);
    sstr = char;
    eidx = strfind(se_in{idx},'e');
    counter = 1;
    while ~isempty(str2num(se_in{idx}(eidx+counter))) || ~isempty(strmatch(se_in{idx}(eidx+counter),'.', 'exact')) || ~isempty(strmatch(se_in{idx}(eidx+counter),'_', 'exact'))
        if se_in{idx}(eidx+counter) == '_'
            estr(counter) = '.';
        else
            estr(counter) = se_in{idx}(eidx+counter);
        end
        if eidx+counter+1 <= size(se_in{idx},2)
            counter = counter+1;
        else
            break
        end
    end
    e(idx) = str2num(estr);
    if e(idx) < s(idx)
        error('all ''e'' numbers must be larger than their ''s'' counterparts.\n you provided s(%d): %f e(%d): %f', idx, s(idx), idx, e(idx));
        return;
    end
    estr = char;
end

% sort by s, then by e
sortmat = sortrows([s' e']);
col1 = sortmat(:,1);
col2 = sortmat(:,2);
s = col1';
e = col2';
    

end