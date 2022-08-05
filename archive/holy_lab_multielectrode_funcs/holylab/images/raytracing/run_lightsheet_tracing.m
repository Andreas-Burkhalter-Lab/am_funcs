optical_params = struct('colR',4.13,...
                        'colCT',3,...
                        'colPosition',5.7,...
                        'colangle',0,...
                        'colmat','sf11',...
                        'cylR',4.13,...
                        'cylCT',3,...
                        'cylPosition',20,...
                        'cylangle',0,...
                        'cylmat','bk7',...
                        'immerse','saline',...
                        'fiberna',0.1,...
                        'nrays',41,...
                        'rayrgb',[1 0 0],...
                        'plotcomponents',1,...
                        'fieldradius',0.5);

fixednames = {'colR','colCT','colPosition'};
% Thorlabs
%fixedvalues = [5.2 2.5 8.4; ...
%               6.2 2.3 10.5];
% Edmund BK7
%fixedvalues = [4.65 1.5 8; ...
%               4.65 2.38 7.41];
% Edmund SF11
fixedvalues = [4.71 2.5 4.6];
nlenses = size(fixedvalues,1);
sigma_lens = zeros(1,nlenses);
for i = 1:nlenses
  for j = 1:length(fixednames)
    optical_params.(fixednames{j}) = fixedvalues(i,j);
  end
  %varynames = {'colPosition','cylPosition'};
  varynames = {'colPosition'};
  lb = [];
  ub = [];
  A = [];
  b = [];
  %A = [1 -1];
  %b = -(optical_params.cylCT + optical_params.colCT);
  %lb = [1 -Inf -Inf];
  %ub = [Inf Inf Inf];
  %fieldnames = {'colPosition'};
  for j = 1:length(varynames)
    x0(j) = optical_params.(varynames{j});
  end

  minops = optimset('Display','iter');
  myfun = @(x) lightsheet_opt(x,optical_params,varynames);
  [x,sigma] = fminsearch(myfun,x0,minops);
  %[x,sigma] = fmincon(myfun,x0,A,b,[],[],lb,ub,[],minops);
  %[x,sigma] = fminunc(myfun,x0,minops);
  sigma_lens(i) = sigma;
  params_best(i,:) = x;
end
return

for i = 1:length(fieldnames)
    optical_params.(fieldnames{i}) = x(i);
end
optical_params.rayrgb = [1 0 0];

figure
subplot(1,2,1)
% Plot components
[c,r] = lightsheet_pcx(optical_params);
for i = 1:length(c); ctmp = c{i}; feval(ctmp{:}); end

% Plot rays & calculate the waist info
[sigma,x2] = lightsheet_traceconfig(optical_params);

% Plot intensity-weighted width vs. x
subplot(1,2,2)
plot(x2,sigma,'k');
ylabel('Std. dev')
xlabel('Position')
