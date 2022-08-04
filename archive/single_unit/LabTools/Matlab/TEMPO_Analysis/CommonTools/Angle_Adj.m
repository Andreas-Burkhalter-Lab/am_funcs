function [pref_adj_one, off_adj] = Angle_Adj(pref_one, pref_two);

if pref_one > 360
    pref_adj_one = mod(pref_one, 360);
elseif pref_one < -360
    pref_adj_one = mod(pref_one, -360);
end

if pref_one < 0
    if pref_one < -180
        pref_adj_one = mod(pref_one, 180);
    else
        pref_adj_one = mod(pref_one, -180);
    end
elseif pref_one > 0
    if pref_one > 180
        pref_adj_one = mod(pref_one, -180);
    else
        pref_adj_one = mod(pref_one, 180);
    end
end

test1 = mod(pref_two, 360);
test2 = mod(pref_two, -360);

if abs(pref_adj_one-test1) < abs(pref_adj_one-test2)
    off_adj = test1;
else
    off_adj = test2;
end