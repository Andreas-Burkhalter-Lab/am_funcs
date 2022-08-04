%% This script will test dynamic range of a certain condition of background noise. 
%% Average activity (density) and the sigma value over time will be recorded for reference.
%% For evoked-->ongoing transients, make r0>r1. For ongoing-->evoked transients, make r0<r1. 
%% For ongoing activity, make r0=r1 and choose small values. For evoked activity, make r0=r1 and choose larger values. 
%% Run for at least 10 cycles in order to eliminate noise. 
%% Synaptic depression parameter u is static. Input goes to all neurons equally. 

tic

run KC_build

cycles = 10

rmin = -5       %% minimum value (exponent of 10) of rstim to test
rint = 0.2      %% intervals of exponents of 10 to test for rstim
rmax = 1.4      %% maximum value (exponent of 10) of rstim to test


r0 = 10^-3           % level of initial low background noise (grey screen)
r1 = 10^-4.5           % level of increased background noise (catcam)
rstim = 10^rmin          % level of input while stimulus (red dot) is applied

    %% Parameters for Synaptic Depression
u0mean = 0.07;
tau = 750;              % tau for W    


T1 = 2000                   % duration of low background noise (grey screen) before background noise increase (catcam)
Tprestim = 10              % duration of increased background noise (catcam) before stimulus is applied
Tstim = 30                  % duration of stimulus - 30 timesteps ~ 90ms from experiment
Tpoststim = 750            % duration of increased background noise (catcam) after stimulus is applied
Ttotal = T1 + Tprestim + Tstim + Tpoststim;

Tresponse = 50;             % time window following the stimulus (red dot) during which activity is measured for stimulus/reponse (dynamic range) 

density = zeros(1,Ttotal);      %% records average activity
sigrec = zeros(1,Ttotal);       %% records sigma during periods of evoked activity

Spikes = zeros(N,1);
X = zeros(N,1);
ii = 0;
c = 0;

for c = 1:cycles
    ii = 0
    for i = rmin:rint:rmax
    ii = ii+1;
    display(i)
    rstim = 10^i;
    r = r0;                  %% set input to the level of ongoing activity
    for t = 1:T1            %% timestep operations 
        Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
        Spikes = mod(Spikes,states-1);
        X = (Spikes==1);                %% logical containing spike/no-spike data
        if ( any(X) ) %& tt<500)       %% if there was a spike during this timestep     
            density(ii,t,c) = sum(Spikes==1)/N;          % record average activity
            sigrec(ii,t,c) = mean(sum(W));                %%% record sigma values

            W = W - u0mean*bsxfun(@times,W,X');             %%% dynamic synapes - short-term depression
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            
         else                                   %% if no spikes fired
            density(ii,t,c) = sum(Spikes==1)/N;          % record average activity
            sigrec(ii,t,c) = mean(sum(W));                %%% record sigma values
            
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
        end
    end

    r = r1;                  %% set input to the level of evoked activity
    for t = T1+1:T1+Tprestim             %% timestep operations      
        Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
        Spikes = mod(Spikes,states-1);
        X = (Spikes==1);                %% logical containing spike/no-spike data
        if ( any(X) ) %& tt<500)       %% if there was a spike during this timestep     
            density(ii,t,c) = sum(Spikes==1)/N;          % record average activity
            sigrec(ii,t,c) = mean(sum(W));                %%% record sigma values

            W = W - u0mean*bsxfun(@times,W,X');             %%% dynamic synapes - short-term depression
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            
         else                                   %% if no spikes fired
            silent = 1;
            density(ii,t,c) = sum(Spikes==1)/N;          % record average activity
            sigrec(ii,t,c) = mean(sum(W));                %%% record sigma values
            
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
        end
    end
	
	r = rstim;               %% set input to stimulus intensity (red dot)
    for t =  T1+Tprestim+1:T1+Tprestim+Tstim            %% timestep operations 
        Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
        Spikes = mod(Spikes,states-1);
        X = (Spikes==1);                %% logical containing spike/no-spike data
        if ( any(X) )                     %% if there was a spike during this timestep     
            density(ii,t,c) = sum(Spikes==1)/N;          % record average activity
            sigrec(ii,t,c) = mean(sum(W));                %%% record sigma values

            W = W - u0mean*bsxfun(@times,W,X');                          %%% dynamic synapes - short-term depression
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            
         else                                   %% if no spikes fired
            silent = 1;
            density(ii,t,c) = sum(Spikes==1)/N;          % record average activity
            sigrec(ii,t,c) = mean(sum(W));                %%% record sigma values
            
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
        end
    end
	
	r = r1;                  %% turn off stimulus and return to background noise 
    for t = T1+Tprestim+Tstim+1:T1+Tprestim+Tstim+Tpoststim           %% timestep operations 
        Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
        Spikes = mod(Spikes,states-1);
        X = (Spikes==1);                %% logical containing spike/no-spike data
        if ( any(X) ) %& tt<500)       %% if there was a spike during this timestep     
            density(ii,t,c) = sum(Spikes==1)/N;          % record average activity
            sigrec(ii,t,c) = mean(sum(W));                %%% record sigma values

            W = W - u0mean*bsxfun(@times,W,X');             %%% dynamic synapes - short-term depression
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            W = W + 0.5 * bsxfun(@times,(W0-W),(~X)') / tau;              %%% dynamic synapses - recovery only for neurons that didn't spike
            
         else                                   %% if no spikes fired
            silent = 1;
            density(ii,t,c) = sum(Spikes==1)/N;          % record average activity
            sigrec(ii,t,c) = mean(sum(W));                %%% record sigma values
            
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
            W = 0.5 * (W0 - W)/tau + W;            %%% dynamic synapses - recovery
        end
    end
end
end

%%%% Dynamic Range Calculations
density = mean(density,3);
sigrec = mean(sigrec,3);

response_window = (density(:,T1+Tprestim+1:Tstim+T1+Tprestim+Tresponse))';
sr = (mean(response_window))';   % stimulus-respones curve

F0 = min(sr);                % uses lowest r and F value for minimum F value
Fmax = max(sr);            % uses final r and F value for Fmax - assumes saturation by final r value
F10 = F0 + 0.1 * (Fmax - F0);
F90 = F0 + 0.9 * (Fmax - F0);
[F10app,F10ind] = min(abs(bsxfun(@minus,sr,F10)))  ;    % F10 = closeness of approximation of F10 to recorded F value, F10ind = index of approximation in sigdens
[F90app,F90ind] = min(abs(bsxfun(@minus,sr,F90))) ;   % F90 = closeness of approximation of F90 to recorded F value, F90ind = index of approximation in sigdens

r10 = 10.^(rint*F10ind + rmin - rint);       % input required to achieve 10% of maximum response
r90 = 10.^(rint*F90ind + rmin - rint);       % input required to achieve 90% of maximum response

delta = 10 * log(r90 ./ r10)                % dynamic range (decibels)

plot(sr)




toc
runtime = toc