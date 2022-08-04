function test2dgauss


rand('state',sum(100*clock));

pars(1) = 20;
pars(2) = 100;
pars(3) = -4;
pars(4) = 5
pars(5) = 0;
pars(6) = 5;

unique_x_pos = [-10 -5 0 5 10];
unique_y_pos = [-10 -5 0 5 10];

%create interpolated arrays for data display
x_interp = unique_x_pos(1): 5 : unique_x_pos(length(unique_x_pos));
y_interp = unique_y_pos(1): 5 : unique_y_pos(length(unique_y_pos));
z_gauss = zeros(length(x_interp), length(y_interp));


num_reps = 3;
raw = zeros(length(x_interp)* length(y_interp)*num_reps,3);
temp = zeros(length(x_interp)* length(y_interp),3);
means = zeros(length(x_interp)* length(y_interp),3);
      

num_trials = length(x_interp) * length(y_interp);

%obtain fitted data for interpolated arrays
for count=1:num_reps
    for i=1:length(x_interp)
        for j = 1:length(y_interp)
            rand('state',sum(100*clock));
            noisy =  gauss2Dfunc(x_interp(i),y_interp(j), pars)* rand;
            cell_index = ((count-1)*(num_trials))+((i-1)*(length(y_interp)) + j);
        
            temp(((i-1)*(length(y_interp)) + j), 1) = x_interp(i);
            temp(((i-1)*(length(y_interp)) + j), 2) = y_interp(j);
            temp(((i-1)*(length(y_interp)) + j), 3) = noisy;
            z_gauss(i,j) = noisy;

        end
        z_gauss;
    end
    raw(cell_index-num_trials+1:cell_index,1) = temp(:,1);
    raw(cell_index-num_trials+1:cell_index,2) = temp(:,2);
    raw(cell_index-num_trials+1:cell_index,3) = raw(cell_index-num_trials+1:cell_index,3) + temp(:,3);
end

means(:,1) = raw(1:num_trials, 1);
means(:,2) = raw(1:num_trials, 2);

add_em_up = zeros(num_trials,1);
for i = 1:num_reps
    add_em_up = raw((i-1)*num_trials+1:(i-1)*num_trials+1+num_trials-1, 3) + add_em_up;
end

means(:,3) = add_em_up/num_reps;


mean_graph = zeros(length(x_interp), length(y_interp));
mean_graph(:, 1) = means(1:5, 3);
mean_graph(:, 2) = means(6:10, 3);
mean_graph(:, 3) = means(11:15, 3);
mean_graph(:, 4) = means(16:20, 3);
mean_graph(:, 5) = means(21:25, 3);

figure
contourf(x_interp, y_interp, mean_graph)
colorbar

mean_graph2 = flipud(mean_graph);
figure
contourf(x_interp, y_interp, mean_graph2)
colorbar

oldpars = pars;
pars = gauss2Dfit(means,raw)
oldpars

x_interp = unique_x_pos(1): 1 : unique_x_pos(length(unique_x_pos));
y_interp = unique_y_pos(1): 1 : unique_y_pos(length(unique_y_pos));
z_gauss = zeros(length(x_interp), length(y_interp));


for i=1:length(x_interp)
    for j = 1:length(y_interp)
        x_interp(i);
        y_interp(j);
        z_gauss(i,j) =  gauss2Dfunc(x_interp(i),y_interp(j), pars);
    end
end

z_gauss = rot90(z_gauss);
figure
contourf(x_interp, y_interp, z_gauss)
colorbar