function [mov0, map, F, meanImg, R0] = generate_shifted_simulation(Npix, ...
    pixpad, ups, Ncells, NT0, radius, sig_noise, xyshift, sigX, sigMIMG, amp0, ...
    neu_amp, mimgSCALE, muf,seedd)
% sigX = 1;
% sigY = 1;

T = size(xyshift,1);
tpad = 5;

ts = 0:1:50;
tau1 = 10;
tau2 = 1;
kernel = exp(-ts/tau1) - exp(-ts/tau2);
kernel = normc(kernel(:));
% plot(kernel)

NT = NT0 + tpad;

% make a background on images, like blood vessels
rng(seedd);

rates = randn(NT + ceil(tau1*10), 10, 10, 'single');
rates = my_conv2(rates, [4 sigX sigX], [1 2 3]);

rates   = max(0, 1 + neu_amp * rates); 
std(rates(:))
rates = bsxfun(@rdivide, rates, mean(rates,1));

mimgneu = my_conv2(randn(10,10), sigMIMG, [1 2]);
% mimgneu = (mimgneu + mimgneu(end:-1:1, end:-1:1) + mimgneu(end:-1:1, :) + mimgneu(:, end:-1:1))/4;

mc      = mean(mean(mimgneu(4:6,4:6)));
mimgneu = mimgneu * sign(mc-mean(mimgneu(:)));
mimgneu = permute(mimgneu, [3 1 2]);

mimgneu = mimgneu - mean(mimgneu(:));
mimgneu = mimgneu/std(mimgneu(:));
mimgneu = single(max(mimgneu+1,0)) + 1;
mimgneu = mean(mimgneu(:)) + mimgSCALE * (mimgneu - mean(mimgneu(:)));

mean(mimgneu(:))

% imagesc(sq(mimgneu))

rates = bsxfun(@times, rates, mimgneu);
%rates = rates - mean(rates(:));
rates = rates / std(rates(:));
% rates = single(max(rates+2, 0)) + 0.5;

R0 = rates;
rates = reshape(filter(kernel, 1, rates(:, :)), size(rates));
rates = permute(rates, [2 3 1]);

mov = single(imresize(rates, [Npix Npix])); 
mov = reshape(mov, Npix^2, []);

icell = 0;

xs = repmat(1:Npix, Npix, 1);
ys = xs';

map = zeros(Npix, 'single');
nCoef = .5 + rand(Ncells,1)/2;
nCoef(:) = .8;

clear xcenter ycenter statGT Fcell
while 1
    % propose new center 
    xcenter(1+icell) = ceil(rand * (Npix-2*radius) + radius);
    ycenter(1+icell) = ceil(rand * (Npix-2*radius) + radius);
    
    ds = (xcenter - xcenter(1+icell)).^2 + (ycenter - ycenter(1+icell)).^2;
    ds = ds(1:icell).^.5;
    
    if icell==0 || min(ds)> radius * 1.33 %2/3
        icell = icell + 1;
        
        dall = (xs - xcenter(icell)).^2 + (ys - ycenter(icell)).^2;
        
        ipix = find(dall.^.5<=radius);
        npix = numel(ipix);
        
        iy = ceil(ycenter(icell) *10/Npix);
        ix = ceil(xcenter(icell) *10/Npix);
        
        R = (rand*2) * mean(R0(:)) * ones(size(R0,1), 1); 
%         R = mean(R0(:)) * ones(size(R0,1), 1); 
                
        sp = muf * poissrnd(R/35);
        
        % to avoid no activity in first bin
        sp(1) = rand;
        
        amp  = amp0 * mean(sq(R0(:, iy, ix)))/mean(R0(:));
        Fcell = amp * filter(kernel, 1, sp(:));
        %keyboard;
        Fcell = bsxfun(@plus,Fcell,mean(Fcell));
        mov(ipix, :) = bsxfun(@plus, nCoef(icell) * mov(ipix, :), Fcell');
                
        map(ipix) =  map(ipix) + 1;
        
%         mimgneu(ipix) = bsxfun(@plus, nCoef(icell) * mimgneu(ipix), mean(Fcell));
        if icell==Ncells
            break;
        end
        
        F(:, icell) = Fcell;
    end
end
mov = mov(:, ceil(tau1*10) + [1:NT]); %

clear xs ys dall;

mov = mov .* max(0, 1+ sig_noise*randn(size(mov), 'single'));
mov = reshape(mov, Npix, Npix, NT);

meanImg = mean(mov,3);
% meanImg = mimgneu;
imagesc(meanImg);
drawnow;

%movShift = zeros(Npix/5,Npix/5,NT,'single');
disp('shifting now');
% shift movie
%keyboard;
mov0 = zeros(Npix/ups,Npix/ups,T,'single');
% loop over movie multiple times with different shifts
for k = 1:(T/NT0)
    for j = 1:NT0
        for d = 1:2
            mshift = round(xyshift(j+(k-1)*NT0,d)*ups);
            indmap{d}  = [max(1,mshift+1) : min(Npix,Npix+mshift)];
            indorig{d} = [max(1,-1*mshift+1) : min(Npix,Npix-mshift)];
        end
        mj = zeros(Npix,Npix,'single');
        mc = mov(:,:,j); % current frame
        mj(indmap{1},indmap{2}) = mc(indorig{1},indorig{2});
        mov0(:,:,j+(k-1)*NT0)             = mj(1:ups:end,1:ups:end);
        if mod(j+(k-1)*NT0,100)==0
            disp(j+(k-1)*NT0);
        end
    end
end
%
mov0        = mov0(pixpad/2+[1:Npix/ups-pixpad],pixpad/2+[1:Npix/ups-pixpad],:);
meanImg     = meanImg(1:ups:end,1:ups:end);
meanImg     = meanImg(pixpad/2+[1:Npix/ups-pixpad],pixpad/2+[1:Npix/ups-pixpad],:);


