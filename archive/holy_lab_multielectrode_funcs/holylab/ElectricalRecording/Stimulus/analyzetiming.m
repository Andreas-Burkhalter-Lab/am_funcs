% This method is crude, but it mostly works for this data set
%mm= [min(swf) max(swf)];
mm = [-0.0149 0.571]; % voltages corresponding to valve 0 and 15, respectively
ssc = round((swf-mm(1))/diff(mm)*15);
%ssc(find(ssc == 1)) = 0;  % There was 1 error made due to voltage fluctuations, fix it
ssc = round(medfilt1(ssc,5,100000));
tswitch = find(ssc(2:end) ~= ssc(1:end-1));
vlv = ssc(tswitch+1);
while (vlv(end) > 0)
  vlv = vlv(1:end-1);
  tswitch = tswitch(1:end-1);
end
indexon = find(vlv);
timeon = (tswitch(indexon+1)-tswitch(indexon))/h.scanrate;
vlvon = vlv(indexon);
[vlvon; timeon]
