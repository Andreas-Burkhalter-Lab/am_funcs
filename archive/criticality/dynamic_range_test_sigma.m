%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%          Dynamic Range (test sigma)
%%%%%%%%%%%%%% This program generates plots of stimulus and response values
%%%%%%%%%%%%%% for different values of sigma to replicate figure 2 of
%%%%%%%%%%%%%% Kinouchi and Copelli (2006).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
tic;
run KC_build
T=1000;
density = zeros(1,T);
Spikes = zeros(N,1);
     
rmin = -5;
rmax = 1;
rint = 0.1;                  % intervals of exponents of 10 for r values
sr = zeros(round((rmax-rmin+rint)/rint),20); % this matrix contains response values for each input and sigma value
hh = 0;

for h=0.1:.1:2
    hh = hh+1;
    sigma=h
%     display(sigma)
    pmax=2*sigma/k;
    v=zeros(1,(N-1)*N/2);                % this vector will be used to form random connections between elements
    v(1,randperm((N-1)*N/2,k*N/2))=1;
    v = pmax * v .* (rand(1,length(v))); % each connection has random strength from 0 to pmax

    W=zeros(N);                          % create weight matrix
    W((triu(ones(N))-eye(N))==1)=v;  
    W=W+W';                              % copy all connections over the diagonal to make them symmetrical
    ii = 0;
    for i = rmin:rint:rmax     % range of exponents of dr to test
        ii = ii+1;
        r=10.^i;
        for t = 2:T
            X = (Spikes==1);
            Spikes = (floor(heaviside((1 - r)*(W * X) + r - rand(N,1)))).* (Spikes==0) + (Spikes~=0).*(1 + Spikes);
            Spikes = mod(Spikes,states-1);
            density (1,t) = sum(Spikes==1)/N;
        end
    sr(ii,hh)= mean(density);
    end 
end
% sr = sr'       ;       % rows are r values, columns are sigma values

F0 = sr(1,:)       ;        % uses lowest r and F value for minimum F value
Fmax = sr(end,:)  ;         % uses final r and F value for Fmax - assumes saturation by final r value
F10 = F0 + 0.1 * (Fmax - F0);
F90 = F0 + 0.9 * (Fmax - F0);
[F10app,F10ind] = min(abs(bsxfun(@minus,sr,F10)))  ;    % F10 = closeness of approximation of F10 to recorded F value, F10ind = index of approximation in sigdens
[F90app,F90ind] = min(abs(bsxfun(@minus,sr,F90))) ;   % F90 = closeness of approximation of F90 to recorded F value, F90ind = index of approximation in sigdens

r10 = 10.^(rint*F10ind - (5 + rint));
r90 = 10.^(rint*F90ind - (5 + rint));

delta= 10 * log(r90 ./ r10);

plot(0.1:0.1:2,delta)
xlabel('sigma')
ylabel('delta (db)')
% 
% plot(-fliplr(-1:rint:5),sigdens)
% % semilogy(-fliplr(-1:rint:5),sigdens)
% xlabel('r(ms^-^1)')
% ylabel('F(ms^-^1)')

runtime = toc;