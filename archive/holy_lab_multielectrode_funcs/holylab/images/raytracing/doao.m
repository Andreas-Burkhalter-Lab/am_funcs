% script to study affect of single deformable mirror 
% adaptive optics
% Diwakar Turaga 2007-01-17

%clear
%close all

% define the rays
r = ray;
r.x0 = [0 0]; % start point
r.e = [1 0]; % start direction
theta = linspace(-pi/8,pi/8,10);
na = 0.5;
nrays = 10;
theta = linspace(-asin(na),asin(na),nrays);
rbase = r;
for i = 1:length(theta);
    r(i) = rbase;
    r(i).e = [cos(theta(i)) sin(theta(i))];
end

% set up of optics
%tissue1
tp1 = struct('normal', [1 0],...
    'c', 0.5, ...
    'apc', [0.5 0],...
    'apr', 5, ...
    'mat1', 'n140',...
    'mat2', 'saline');

%tissue2
tp2 = struct('normal', [1 0],...
    'c', 1, ...
    'apc', [1 0],...
    'apr', 5, ...
    'mat1', 'saline',...
    'mat2', 'air');

%perfect optics
pp = struct('xc', [5 0], ...
    'normal', [1 0], ...
    'f', [3 3], ...
    'plotsize', 6);

%dm1 = struct('mirEdge', [5 -6.5; 16 4.5]);
dm1 = struct('mirEdge', [6 -6; 6 6], ...
             'ml', 0.5);

%DM2
%dm2 = struct('mirEdge', [5 -6; 16 5]);
dm2 = struct('mirEdge', [5.5 -6; 5.5 6]);

%screen
s = struct('normal', [1 0], ...
    'c', 20, ...
    'apc', [20 0], ...
    'apr', 6, ...
    'mat1', 'air', ...
    'mat2', 'air');

%trace optics
figure;
c{1} = opt2dline(tp1);
c{2} = opt2dline(tp2);
c{3} = opt2dperfect(pp);
c{4} = opt2dDM(dm1);
c{5} = opt2dDM(dm2);
c{6} = opt2dline(s);
axis equal;

%trace rays
c{1} = {@opt2dline,tp1};
c{2} = {@opt2dline,tp2};
c{3} = {@opt2dperfect, pp};
c{4} = {@opt2dDM, dm1};
c{5} = {@opt2dDM,dm2};
c{6} = {@opt2dline,s};

for i = 1:length(r)
     routput(i) = raytrace(r(i),c,[1 0 0]); % see what comes out before hitting the screen
end


%trace through the first DM to find point of intersection
for i = 1:length(r)
    rf(i) = raytrace(r(i),c(1:4),[1 0 0]);
end

[err,x0] = ray_convergence_error(rf);

% now trace through just before DM and use the point to deform the mirror
for i = 1:length(r)
    rf(i) = raytrace(r(i),c(1:3),[1 0 0]);
end  

dm1.p = x0;

for rc = 1:length(rf)
    [dmnew] = DMmanipulator(dm1,rf(rc));
    dm1new(2*rc-1,:) = dmnew(1,:);
    dm1new(2*rc,:) = dmnew(2,:);
end

% found the deformable mirror components

% now actual test
dm1.mirEdge = dm1new;

figure;        
opt2dline(tp1);
opt2dline(tp2);
opt2dperfect(pp);
opt2dDM(dm1);
opt2dDM(dm2);
opt2dline(s);
axis equal;

%trace rays
c{1} = {@opt2dline,tp1};
c{2} = {@opt2dline,tp2};
c{3} = {@opt2dperfect, pp};
c{4} = {@opt2dDM, dm1};
c{5} = {@opt2dDM,dm2};
c{6} = {@opt2dline,s};

%all the above to find perfect DM 

% now scan through various points in tissue to find the beam waist for each
% of the points
tisDim = [0 -.2; 0.5 .2];
xsize = 0.01; % step size
ysize = 0.01; % step size
count = 1;
for t1 = tisDim(1,1):xsize:tisDim(2,1)
    for t2 = tisDim(1,2):ysize:tisDim(2,2)

        r = ray;
        r.x0 = [t1 t2]; % start point
        %r.x0 = [0 0];
        r.e = [1 0]; % start direction
        %theta = linspace(-pi/8,pi/8,10);
        theta = linspace(-asin(na),asin(na),nrays);
        rbase = r;
        for i = 1:length(theta);
            r(i) = rbase;
            r(i).e = [cos(theta(i)) sin(theta(i))];
        end


        count2 = 1;
        for i = 1:length(r)
            rou = raytrace(r(i),c,[1 0 0]); % see what comes out before hitting the screen
            %rout = raytrace(r(i),c);
            if isfield(rou,'x0')
                rout(count2) = rou;
                count2 = count2+1;
            end
        end

        [err,x0] = ray_convergence_error(rout);
        
        min_x = min(x0);
        min_waist = min(err);
        
        min_values(count,:) =[t1 t2 min_x min_waist];
        
        count=count+1;
        
        %close all
    end
end

count = 1;
c1 = 1;
for t1 = tisDim(1,1):xsize:tisDim(2,1)
    c2 = 1;
    for t2 = tisDim(1,2):ysize:tisDim(2,2)
        min_waist_grid(c2, c1) = min_values(count,4);
        count = count+1;
        c2 = c2+1;
    end
    c1 = c1+1;
end

xGridRange = [tisDim(1,1) tisDim(2,1)];
yGridRange = [tisDim(1,2) tisDim(2,2)];

figure;
% 
% subplot(1,2,2);
% plot(min_values(:,1), min_values(:,4),'o');
% xlabel('depth (mm)');
% ylabel('error');
% title('error vs. depth')
% 
% subplot(1,2,1);

imagesc(xGridRange,yGridRange,sqrt(min_waist_grid));
xlabel('depth (mm)');
ylabel('length (mm)');
title('length vs. depth with error');
%axis equal
