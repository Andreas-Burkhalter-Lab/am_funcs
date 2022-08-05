function [U,M] = maxproj(Z)
  [Q,R] = qr(Z,0);
  theta = fminsearch(@(th) maxproj_proj(th,Q,Z),0)
  rot = [cos(theta) sin(theta); -sin(theta) cos(theta)];
  U = Q*rot';
  M = rot*R;
  theta = linspace(0,2*pi,100);
  for i = 1:length(theta)
    v(i) = -maxproj_proj(theta(i),Q,Z);
  end
  plot(theta,v)
  
function p = maxproj_proj(theta,Q,Z)
  rot = [cos(theta) sin(theta); -sin(theta) cos(theta)];
  U = Q*rot';
  p = -sum(sum(U.*Z).^2);
  