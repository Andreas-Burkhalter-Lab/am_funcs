%[c,r] =
%objective_aspheric(struct('wd',2.25,'immerse','air','lambda',510));
[c,r] = objective_aspheric(struct('wd',2.5,'immerse','water','lambda',510,'na',0.44));
cla
for i = 1:length(c)
  ctmp = c{i};
  feval(ctmp{:});
end

clear rf
for i = 1:length(r)
  rf(i) = raytrace(r(i),c,[1 0 0]);
end
% Trace the rays to intersection with optic axis
clear xf
for i = 1:length(r)
   xf(i) = rf(i).x0(1) - rf(i).x0(2)*rf(i).e(1)/rf(i).e(2);
   line([rf(i).x0(1) xf(i)],[rf(i).x0(2) 0],'Color',[0 1 0]);
end
sinthetamax = rf(end).e(2)
errors = (max(xf)-min(xf))*sinthetamax


axis equal
shg