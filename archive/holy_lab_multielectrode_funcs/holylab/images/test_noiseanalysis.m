% A test script for imnoise analysis
% Note a tricky issue in imnoiseanalysis: in simple=0 cases, noise is
% gaussian. In simple=1 cases, noise is poisson. Pay attention to how this
% script simulates noise, or you may be surprised by discrepancies.

G = [1];
t = [0.1 0.2 0.5 1 2 5];
%G = 1;
%t = 1;
nreps = 50;
nx = 40;
ny = 30;
% The signal
rsigt = rand(ny,nx)*100;
nf = 3;
bf = 5;
toclip = 1;
imrange = [0 400];

% ngf = 0;
% bgf = 0;
% ncf = 0;
% bcf = 0;

% Gain settings
gaint = [];
for i = 1:length(G);
  %gaint(i) = rand(1) + G(i);
  gaint(i) = G(i);
end
%biast = rand(ny,nx)*bf;
biast = bf;
noiset = nf;
ip = repmat(struct,[1 length(t)*length(G)*nreps]);
counter = 1;
for i = 1:length(t);
  for j=1:length(G);
    for k=1:nreps;
      ip(counter).image = gaint(j).*(poissrnd(rsigt*t(i))) + biast + randn(ny,nx)*noiset;
      if (toclip)
        indxbad = find(ip(counter).image > imrange(end));
        ip(counter).image(indxbad) = imrange(end);
      end
      ip(counter).gainsetting = G(j);
      ip(counter).inttime = t(i);
      ip(counter).imrange = imrange;
      counter = counter+1;
      if (mod(counter,100) == 0)
        counter
      end
    end;
  end;
end
for k = 1:length(ip)
  ip(k).xrange = [1 nx];
  ip(k).yrange = [1 ny];
end
% [rsigm,noisegm,biasgm,noisecm,biascm,gainm] =
% imnoiseanalysis(ip,struct('plotraw',1));