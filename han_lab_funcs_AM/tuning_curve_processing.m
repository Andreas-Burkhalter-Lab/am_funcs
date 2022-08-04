function f = tuning_curve_processing(feature_count, prior_prob)
winsize=64;
alpha=50;
gkern=gausswin(winsize,alpha);
%Also try with gaussian filter, see what works better
gf=fspecial('gaussian', [32 1], .5);

% feature_count = bsxfun(@rdivide, feature_count, prior_prob');

% feature_countgf = conv2(feature_count, gkern/sum(gkern),'same'); %Widest
feature_count=reshape(feature_count(~isnan(feature_count)),[],size(feature_count,2));
feature_countgf = imfilter(feature_count, gf); %Sharpest

% feature_countnc = normc(feature_count);
feature_countn = bsxfun(@rdivide, feature_countgf, max(feature_countgf));

% feature_count=feature_countn;

f=feature_countn;

% figure()
% hold on
% plot(feature_count(:,50),'ro')
% plot(feature_countgc(:,50),'g')
% plot(feature_countgf(:,50),'y')


% figure()
% imagesc(feature_count)
% title('smoothed')
% colorbar()
% figure()
% imagesc(feature_countn)
% title('smoothed and normed')
% colorbar()
% 
% figure()
% imagesc(feature_countnc)
% title('smoothed and normced')
% colorbar()
% return f
end