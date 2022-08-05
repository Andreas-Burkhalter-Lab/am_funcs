function filter = punctate_filter_jpm_1()
% this filter is in process, just trying to find a filter which will enhance bright cells of the size seen in FJC expts.

base_unit = 0.027777; % since values need to = 1 to keep relative values constant and values add to 96*base_unit
bu = base_unit; % for shorthand
filter = zeros(14);
filter(1,:)  = -bu;
filter(2,:)  = -bu;
filter(3,:)  = [-bu   -bu   -bu   0     0     0     0     0     0     0     0     -bu   -bu   -bu];
filter(4,:)  = [-bu   -bu   0     0     0     2*bu  2*bu  2*bu  2*bu  0     0     0     -bu   -bu];
filter(5,:)  = [-bu   -bu   0     0     2*bu  2*bu  3*bu  3*bu  2*bu  2*bu  0     0     -bu   -bu];
filter(6,:)  = [-bu   -bu   0     2*bu  2*bu  3*bu  3*bu  3*bu  3*bu  2*bu  2*bu  0     -bu   -bu];
filter(7,:)  = [-bu   -bu   0     2*bu  3*bu  3*bu  4*bu  4*bu  3*bu  3*bu  2*bu  0     -bu   -bu];
filter(8,:)  = [-bu   -bu   0     2*bu  3*bu  3*bu  4*bu  4*bu  3*bu  3*bu  2*bu  0     -bu   -bu];
filter(9,:)  = [-bu   -bu   0     2*bu  2*bu  3*bu  3*bu  3*bu  3*bu  2*bu  2*bu  0     -bu   -bu];
filter(10,:) = [-bu   -bu   0     0     2*bu  2*bu  3*bu  3*bu  2*bu  2*bu  0     0     -bu   -bu];
filter(11,:) = [-bu   -bu   0     0     0     2*bu  2*bu  2*bu  2*bu  0     0     0     -bu   -bu];
filter(12,:) = [-bu   -bu   -bu   0     0     0     0     0     0     0     0     -bu   -bu   -bu];
filter(13,:) = -bu;
filter(14,:) = -bu;

% sum(sum(filter,2),1)  % for testing

% figure(2); imshow(filter); set(get(2, 'children'), 'clim', [-base_unit 4*base_unit]);  % for testing 