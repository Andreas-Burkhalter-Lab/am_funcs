 %%%% This program tests the duration of one avalanche, starting with a
 %%%% fraction of the neurons active, for different network sizes. Synapses
 %%%% are static.

tic;
run KC_build;

init = 0.1;         % fraction of neurons starting in the active state
Spikes = floor(heaviside(rand(N,1)-(1-init)));            % vector representing the activity of all neurons at a time; start with 10% activity
density = zeros(1,T);            % density of active sites at each timepoint
density(1,1) = sum(Spikes)/N;
 
 duration = zeros(100,100);
N=0

for i = 1:100   
    for j=1:100     % j = number of samples to average
    N = 30*i;
  
    v=zeros(1,(N-1)*N/2);                % this vector will be used to form random connections between elements
    v(1,randperm((N-1)*N/2,k*N/2))=1;
    v = pmax * v .* (rand(1,length(v))); % each connection has random strength from 0 to pmax
    W=zeros(N);                          % create weight matrix
    W((triu(ones(N))-eye(N))==1)=v;  
    W=W+W';                              % copy all connections over the diagonal to make them symmetrical
   
    Spikes = zeros(N,1);
    tt=0;
    Spikes = floor(heaviside(rand(N,1)-(1-init)));
    while (any(Spikes==1) )
        X = (Spikes==1);
        Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
        %        Spikes = (floor(heaviside(W * X -rand(N,1) )) | (floor(heaviside(input(:,t) - rand(N,1))))) .* (Spikes==0) + (Spikes~=0).*(1 + Spikes); % see 'input' above                                                                                                                                                              
        Spikes = mod(Spikes,states-1);
        tt = tt+1;
    end
    duration(j,i) = tt;
    end
end

plot(30:30:3000, mean(duration))
xlabel('Number of Elements N')
ylabel('Mean Duration (ms)')