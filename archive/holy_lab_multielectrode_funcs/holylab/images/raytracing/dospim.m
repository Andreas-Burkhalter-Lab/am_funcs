% Build the components: make sure that the one you want is the last
% uncommented one
% Do the whole thing with plano-convex lenses
optical_params = struct('colR',5.2,...
                        'colCT',2.5,...
                        'colPosition',[6.4 0],...
                        'colangle',0,...
                        'colmat','bk7',...
                        'cylR',4.13,...
                        'cylCT',3,...
                        'cylPosition',[20 0],...
                        'cylangle',0,...
                        'cylmat','bk7',...
                        'immerse','saline',...
                        'fiberna',0.1,...
                        'nrays',41,...
                        'rayrgb',[1 0 0],...
                        'fieldradius',0.5);

clf
subplot(1,2,1)

% Plot components
[c,r] = lightsheet_pcx(optical_params);
for i = 1:length(c); ctmp = c{i}; feval(ctmp{:}); end

% Plot rays & calculate the waist info
[sigma,x] = lightsheet_traceconfig(optical_params);

% Plot intensity-weighted width vs. x
subplot(1,2,2)
plot(x,sigma,'k');
ylabel('Std. dev')
xlabel('Position')
