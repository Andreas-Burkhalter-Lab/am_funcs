[c,r] = objective_grin(struct('wd','infinity-focused','immerse','saline','y',0,'lambda', ...
			      510));
traceback = false;

cla
for i = 1:length(c)
  ctmp = c{i};
  feval(ctmp{:});
end

clear rf rf_match;
rfIndex = 0;
for i = 1:length(r)
  rtmp = raytrace(r(i),c,[1 0 0]);
  if (length(fieldnames(rtmp)) > 0)
    rf(rfIndex+1) = rtmp;
    rf_match(rfIndex+1) = i;
    rfIndex = rfIndex+1;
  end
end

clear xf
if (length(c) < 3)
  % Trace the rays coming from GRIN to their intersection with the optic axis
  for i = 1:length(rf_match)
    %i = rf_match(ii);
    xf(i) = rf(i).x0(1) - rf(i).x0(2)*rf(i).e(1)/rf(i).e(2);
    if traceback
      line([rf(i).x0(1) xf(i)],[rf(i).x0(2) 0],'Color',[0 1 0]);
    end
  end
  % Calculate NAout
  NAout = abs(rf(1).e(2))
  % Calculate the significance of the z-errors
  sinthetamax = rf(end).e(2)
  zerrors = (max(xf)-min(xf))
  proderrors = zerrors*sinthetamax
end

% Calculate NAin
e = r(rf_match(1)).e;
NAin = opt_refrindx(c{1}{2}.mat1,r(1).w) * abs(e(2))

axis equal
shg