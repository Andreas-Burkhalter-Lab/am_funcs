%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%          Dynamic Range - Test Background Noise
%%%%%%%%%%%%%% This program generates plots of stimulus-functions for different
%%%%%%%%%%%%%% levels of background noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
tic
run KC_build

T = 1000;
density = zeros(1,T);
rmin = -5;
rmax = 0.2;
rint = .1;
sr = zeros(11,round((rmax-rmin+rint)/rint));
Spikes = zeros(N,1);
r1 = 0;                 %% background noise

jj = 0;
    for j=-1:.1:0            % values of r0 background noise to test
        jj = jj+1;
        ii = 0;
        r0 = 10.^j;
%         display(j)
        for i = rmin:rint:rmax     % range of exponents of dr to test
            ii = ii+1;
            r1=10.^i;
            r = r0 + r1;
            for t = 1:T
                X = (Spikes==1);
                Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
                Spikes = mod(Spikes,states-1);
                density (1,t) = sum(Spikes==1)/N;
            end
            sr(ii,jj) = mean(density);
        end 
    end

F0 = sr(1,:)       ;        % uses lowest r and F value for minimum F value
Fmax = sr(end,:)  ;         % uses final r and F value for Fmax - assumes saturation by final r value
F10 = F0 + 0.1 * (Fmax - F0);
F90 = F0 + 0.9 * (Fmax - F0);
[F10app,F10ind] = min(abs(bsxfun(@minus,sr,F10)))  ;    % F10 = closeness of approximation of F10 to recorded F value, F10ind = index of approximation in sigdens
[F90app,F90ind] = min(abs(bsxfun(@minus,sr,F90))) ;   % F90 = closeness of approximation of F90 to recorded F value, F90ind = index of approximation in sigdens

r10 = 10.^(rint*F10ind - (5 + rint));
r90 = 10.^(rint*F90ind - (5 + rint));

delta= 10 * log(r90 ./ r10);

plot(rmin:rint:rmax,sr)

toc