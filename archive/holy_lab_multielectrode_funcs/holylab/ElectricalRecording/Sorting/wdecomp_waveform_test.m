n_samples = 1000000;
n_events = 50000;
w_true = [0 1 2 1.5 1 0.5 0 -0.5 -1 0; 3 2 3 2 3 2 0 -2 2 3];
%w_true = [0 1 2 1.5 1 0.5 0 -0.5 -1 0];
noise_amplitude = 0.1;
trange = [-2 7];

n_components = size(w_true,1);
a = randn(n_components,n_events);
%a = ones(1,n_events);  % uncomment if you want uniform event sizes

t = randperm(n_samples);
t = sort(t(1:n_events));
keepFlag = t + trange(1) > 0 & t + trange(2) <= n_samples;
t = t(keepFlag);
a = a(:,keepFlag);
n_w = diff(trange)+1;

%t = n_w:2*n_w:n_samples-n_w;  % uncomment this to insure no overlaps

v = zeros(1,n_samples);
rng = trange(1):trange(2);
for i = 1:length(t)
  v(t(i)+rng) = v(t(i)+rng) + a(:,i)'*w_true;
end
v = v + noise_amplitude * randn(1,n_samples);

%w_test = wdecomp_waveform_singlecomp(v,t,trange,a);
w_test = wdecomp_waveform(v,t,trange,a);
w_test = squeeze(w_test)';
