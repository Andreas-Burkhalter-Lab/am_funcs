function vol = zern2mirao1(n,m,w, dm_angle, max_vol)

% zern2mirao1 is the dumb way to optimize actuator voltages
% it will not take into consideration influence functions


% input will be 
% n: order N of zernike
% m: frequency M
% w: weight of the particular zernike term
% dm_angle: tilt of the DM, currently pi/4 for angle between incident beam 
% and mirror orthogonal axis

% output will be a vector of size 52, containing voltages for each of the
% 52 mirror lets

% n = [0  1  1  2  2  2  3  3  3  3];
% m = [0 -1  1 -2  0  2 -3 -1  1  3];
% w = [0 0 0 0 1 0 0 0 0 0];
% dm_angle = pi/4;
% max_vol = 0.5

% Diwakar Turaga May 1st 2008

act_loc = [NaN NaN 11 19 27 35 NaN NaN; ...
           NaN 5 12 20 28 36 43 NaN; ...
           1 6 13 21 29 37 44 49; ...
           2 7 14 22 30 38 45 50; ...
           3 8 15 23 31 39 46 51;...
           4 9 16 24 32 40 47 52;...
           NaN 10 17 25 33 41 48 NaN; ...
           NaN NaN 18 26 34 42 NaN NaN];

% there are two ways to do this
% first is the dumbway
% calculate only at the actuator position
% no influence function stuff is used

dm_grid_x = -8.75:2.5:8.75; % in mm
dm_grid_y = dm_grid_x*sin(dm_angle);

[X,Y] = meshgrid(dm_grid_x,dm_grid_y);
[theta,r] = cart2pol(X,Y);

r = r/max(max(r)); %normalizing the size of the circle


    
idx = r<=1;
z = nan(size(X));

y = zernfun(n,m,r(idx),theta(idx), 'norm');

z(idx) = y(:,1);
z_final = z;

for i = 2:length(n)
    z(idx) = y(:,i);
    z_final = z_final+ z*w(i);
end

z_final = z_final*max_vol/max(max(abs(z_final)));

for j = 1:52
    act_loc_j = act_loc == j;
    act_vol_grid = act_loc_j.*z_final;
    vol(j) = sum(sum(act_vol_grid));
end

figure; imagesc(z_final)
set(gca,'XTick',[],'YTick',[])
axis image; colorbar





















    

    

