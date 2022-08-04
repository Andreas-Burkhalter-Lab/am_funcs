%%% Use this program to create data for avalanches during stable ongoing,
%%% stable evoked, evoked-->ongoing transient, and ongoing--evoked
%%% transient activity. See the section labeled 'TRANSIENT AVALANCHES' for
%%% plotting pdfs of transient and non-transient avalanches. 

tic

run KC_build


r0 = 1 * 10^-4.5;       %% level of ongoing activity input
r1 = 1 * 10^-3;         %% level of evoked activity input


    %% Parameters for Synaptic Depression
tau = 750;              % tau for W    
u0mean = 0.07;

Spikes = zeros(N,1);
X = zeros(N,1);

fired = Spikes;
firings = [];
cycles = 100;                 %% number of times to run ongoing and evoked activity
T1 = 3000;                   %% number of timesteps spent in ongoing activity per cycle
T2 = 3000;                   %% number of timesteps spent in evoked activity per cycle

Ttotal = cycles * (T1 + T2);
Burst_size = zeros(Ttotal,4);

densong0 = zeros(1,T1);          %% records average activity during the period of ongoing activity
densevo = zeros(cycles,T2);      %% records average activity during periods of evoked activity
densong1 = zeros(cycles,T1);    %% records average activity during periods of ongoing activity after the first

sigong0 = zeros(1,T1);          %% records sigma during the period of ongoing activity
sigevo = zeros(cycles,T2);       %% records sigma during periods of evoked activity
sigong1 = zeros(cycles,T1);     %% records sigma during periods of ongoing activity after the first

silent = 1;
avtime = 0;
avcycle = 1;
av_evoked = 0;

for c = 1:cycles
    display(c)
    r = r0;                  %% set input to the level of ongoing activity
    tt = 0;
    for t = 1:T1            %% timestep operations 
        tt = tt+1;           %% counts absolute number of timesteps in this simulation    
        Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
        Spikes = mod(Spikes,states-1);
        X = (Spikes==1);                %% logical containing spike/no-spike data
        if ( any(X) )                  %% if there was a spike during this timestep     

            if (silent == 1)               %% detects whether or not this timestep is the first timestep of an avalanche
                av_evoked = 0;
                avcycle = c;
                avtime = t;             %% records the time of the beginning of this avalanche to be added to the second column of Burst_size
                silent = 0 ;             %% set 'silent' off so as not to trigger again at the next timestep of the avalanche
            end
            densong1(c,t) = sum(Spikes==1)/N;          % record average activity
            sigong1(c,t) = mean(sum(W));                %%% record sigma values
            W = W - u0mean*bsxfun(@times,W,X');             %%% dynamic synapes - short-term depression
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            
            fired = fired | X;
         else                                   %% if no spikes fired
            silent = 1;
            densong1(c,t) = sum(Spikes==1)/N;          % record average activity
            sigong1(c,t) = mean(sum(W));                %%% record sigma values
            
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
            
            Burst_size(tt,:) = [sum(fired),av_evoked,avcycle,avtime];
            fired = zeros(N,1);
        end
    end
    
    r = r1;                  %% set input to the level of evoked activity
    for t = 1:T2            %% timestep operations 
         tt = tt+1;           %% counts absolute number of timesteps in this simulation
        Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
        Spikes = mod(Spikes,states-1);
        X = (Spikes==1);                %% logical containing spike/no-spike data
        if ( any(X) )                  %% if there was a spike during this timestep     

            if (silent == 1)               %% detects whether or not this timestep is the first timestep of an avalanche
                av_evoked = 1;
                avcycle = c;
                avtime = t;             %% records the time of the beginning of this avalanche to be added to the second column of Burst_size
                silent = 0 ;             %% set 'silent' off so as not to trigger again at the next timestep of the avalanche
            end
            densevo(c,t) = sum(Spikes==1)/N;          % record average activity
            sigevo(c,t) = mean(sum(W));                %%% record sigma values

            W = W - u0mean*bsxfun(@times,W,X');             %%% dynamic synapes - short-term depression
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            
            fired = fired | X;
         else                                   %% if no spikes fired
            silent = 1;
            densevo(c,t) = sum(Spikes==1)/N;          % record average activity
            sigevo(c,t) = mean(sum(W));                %%% record sigma values
            
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery

            Burst_size(tt,:) = [sum(fired),av_evoked,avcycle,avtime];
            fired = zeros(N,1);
        end
    end
end

densong0 = densong1(1,:);
densong1 = densong1(2:end,:);
sigong0 = sigong1(1,:);
sigong1 = sigong1(2:end,:);

mean_densevo = mean(densevo,1);
mean_densong1 = mean(densong1,1);
mean_sigevo = mean(sigevo,1);
mean_sigong1 = mean(sigong1,1);

burstevo(:,4) = xburstevo(:,4) + T1;
burstong1(:,4) = xburstong1(:,4) + T1+T2;

tc_dens = [densong0, mean_densevo, mean_densong1];
tc_sig = [sigong0,  mean_sigevo, mean_sigong1];

%%% The following lines generates a scatter plot of avalanches from all three periods. 
figure
scatter([burstong0(:,4) ; burstevo(:,4) ; burstong1(:,4)],[burstong0(:,1) ; burstevo(:,1) ; burstong1(:,1)]);
xlabel('time');
ylabel('Burst Size')

toc
runtime = toc


            %%% TRANSIENT AVALANCHES %%%
%%% Use this section to generate a pdf for specific periods within evoked or ongoing activity 
%%% Adjust the variable 'transient' to decide the duration of time in which
%%% to consider transients after a change in background noise. 
transient = 500; 

%% Plot transient or non-transient avalanches during evoked activity.
tburstevo = xburstevo( (xburstevo(:,4) < (transient+T1)) , :);   %% To plot avalanches during transients, use < here. To plot avalanches outside of transients, use >.
tL_evo = max(tburstevo(:,1));
ts_evo = 1:tL_evo;
[tpdf_evo junk] = histc(tburstevo(:,1),ts_evo);

figure
loglog(ts_evo,tpdf_evo/sum(tpdf_evo),'r','LineWidth',1.1);
xlabel('Size', 'Fontsize',15);
ylabel('PDF', 'Fontsize',15);
title('Probability Size Distribution - Evoked', 'Fontsize',17); 
 

%% Plot transient or non-transient avalanches during ongoing activity. 
tburstong1 = xburstong1( (xburstong1(:,4) < (transient+T1+T2)) , :);   %% To plot avalanches during transients, use < here. To plot avalanches outside of transients, use >.
tL_ong1 = max(tburstong1(:,1));
ts_ong1 = 1:tL_ong1;
[tpdf_ong1 junk] = histc(tburstong1(:,1),ts_ong1);

figure
loglog(ts_ong1,tpdf_ong1/sum(tpdf_ong1),'r','LineWidth',1.1);
xlabel('Size', 'Fontsize',15);
ylabel('PDF', 'Fontsize',15);
title('Probability Size Distribution - Ongoing', 'Fontsize',17);  

