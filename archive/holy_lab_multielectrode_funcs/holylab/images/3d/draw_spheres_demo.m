roic = [0 0 0; 3.4 0 0; 2 5 1; 0 -2 0]; % ROI centers
clust = [1 1 2 3];  % cluster labels
col = [1 0 0; 0 1 0; 0 0 1];  % cluster colors
r = 0.7;  % sphere radius
save_movie = false;

% Draw the spheres
figure; hold on
[xs,ys,zs] = sphere(20);
xs = r*xs; ys = r*ys; zs = r*zs;
n_rois = size(roic,1);
hs = zeros(1,n_rois);
for i = 1:n_rois
  hs(i) = surf(xs+roic(i,1),ys+roic(i,2),zs+roic(i,3));
  set(hs(i),'FaceColor',col(clust(i),:));
end
% Set transparency, lighting
set(hs,'EdgeColor','none',...
  'FaceAlpha',0.6,...
  'FaceLighting','phong',...
  'BackFaceLighting','unlit',...
  'AmbientStrength',0.5)
axis equal
axis vis3d off
lighting phong
light % can specify position, e.g., light('Position',[0 -2 1]), or make lights from multiple directions

% Now play some fun games with view angle---see also campos, camva, camup,
% etc.
frame_counter = 0;
campos([0 0 75])
if save_movie, fn = sprintf('movie%d',frame_counter); frame_counter=frame_counter+1; print('-dpng',fn); else, pause(0.05); end
dth = 3; % angle change, in degrees
for i = 1:90/dth
  camorbit(0,dth);
  if save_movie, fn = sprintf('movie%d',frame_counter); frame_counter=frame_counter+1; print('-dpng',fn); else, pause(0.05); end
end
for i = 1:45/dth
  camorbit(dth,0)
  if save_movie, fn = sprintf('movie%d',frame_counter); frame_counter=frame_counter+1; print('-dpng',fn); else, pause(0.05); end
end
for i = 1:90/dth
  camorbit(-dth,0)
  if save_movie, fn = sprintf('movie%d',frame_counter); frame_counter=frame_counter+1; print('-dpng',fn); else, pause(0.05); end
end

% Assemble into a movie with this command (from the command line):
%  ffmpeg -r 10 -i movie*.png -qscale 5 -y -an test.avi
% and then watch with
%  ffplay test.avi
