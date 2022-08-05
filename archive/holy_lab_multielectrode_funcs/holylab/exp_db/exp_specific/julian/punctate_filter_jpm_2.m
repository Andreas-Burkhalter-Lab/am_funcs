function filter = punctate_filter_jpm_2()
% this filter is in process, just trying to find a filter which will enhance bright cells of the size seen in FJC expts.

base_unit = 0.08333; % since values need to = 1 to keep relative values constant and values add to 96*base_unit
bu = base_unit; % for shorthand
filter = zeros(14);
filter(1,:)  = [-bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu];
filter(2,:)  = [-bu   -bu   -bu   -bu   -bu   1*bu  1*bu  1*bu  1*bu  -bu   -bu   -bu   -bu   -bu];
filter(3,:)  = [-bu   -bu   -bu   -bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu   -bu   -bu   -bu];
filter(4,:)  = [-bu   -bu   -bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu   -bu   -bu];
filter(5,:)  = [-bu   -bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu   -bu];
filter(6,:)  = [-bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu];
filter(7,:)  = [-bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu];
filter(8,:)  = [-bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu];
filter(9,:)  = [-bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu];
filter(10,:) = [-bu   -bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu   -bu];
filter(11,:) = [-bu   -bu   -bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu   -bu   -bu];
filter(12,:) = [-bu   -bu   -bu   -bu   1*bu  1*bu  1*bu  1*bu  1*bu  1*bu  -bu   -bu   -bu   -bu];
filter(13,:) = [-bu   -bu   -bu   -bu   -bu   1*bu  1*bu  1*bu  1*bu  -bu   -bu   -bu   -bu   -bu];
filter(14,:) = [-bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu   -bu];

sum(sum(filter,2),1)  % for testing

%figure(2); imshow(filter); set(get(2, 'children'), 'clim', [-bu 1*bu]);