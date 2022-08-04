function [theta_max, theta_min, max_depth, min_depth] = grad_angle(ap_size, m_disp, grad_mag)

view_dist = 57;
inter_oc = 3;
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};

for i=1:length(ap_size)
   figure
   for j=1:length(grad_mag)
      delta_disp = grad_mag(j) * .5 * ap_size(i);
      max_disp = m_disp + delta_disp;
      min_disp = m_disp - delta_disp;

      max_depth = depth_from_disp(max_disp, view_dist, inter_oc);
      min_depth = depth_from_disp(min_disp, view_dist, inter_oc);
      %max_depth is relative to viewing distance
      max_depth = max_depth - view_dist
      min_depth = min_depth - view_dist;

      distance_max = sqrt(power(ap_size/2,2) + power(max_depth,2))
      distance_min = sqrt(power(ap_size/2,2) + power(min_depth,2));

      theta_max = (acos((ap_size/2)/distance_max))
      theta_max = theta_max * 180/pi
      theta_min = (acos((ap_size/2)/distance_min));
      theta_min = theta_min * 180/pi;

      hold on
      plot([max_depth min_depth], [ap_size(i)/2 -ap_size(i)/2], lines{j});
      grid on
      s{j} = sprintf('%3.2f', theta_max);
   end
   legend(s, 0)
   hold on
   rectangle('Position', [-(ap_size/4) -(ap_size/4) ap_size/2 ap_size/2], 'Curvature', [1 1])
   axis equal
end