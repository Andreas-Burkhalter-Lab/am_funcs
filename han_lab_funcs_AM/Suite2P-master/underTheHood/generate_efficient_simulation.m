function [mALL, map, ipix, F, spall, meanImg, R0, xyshift] = generate_efficient_simulation(Npix,  pixpad,...
    Ncells, NT, radius, sig_noise, sigX, sigMIMG, amp0, ...
    neu_amp, mimgSCALE, muf)

% This simulation is complex, and not commented. 
seedd           = [999 10232 1231 8472 1239 9999 121 230 777 127712 4492 823];
seedd = seedd(2);

rng(5);
xyshift       = 20*exprnd(.4,2,NT).*(rand(2,NT)-.5);%.*(rand(2,NT)<.15);
xyshift       = my_conv2(xyshift',2, 1);
xyshift       = bsxfun(@minus,xyshift,mean(xyshift,1));
xyshift       = round(xyshift);

ts = 0:1:50;
tau1 = 10;
kernel = exp(-ts/tau1);
kernel = normc(kernel(:));

ntpad = 100;

% make a background on images, like blood vessels
rng(seedd);

rates = randn(NT+ntpad, 10, 10, 'single');
rates = my_conv2(rates, [10 sigX sigX], [1 2 3]);

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
R0 = rates;
R0 = R0/std(R0(:));
R0 = max(0, R0 - 1);

rates = bsxfun(@times, rates, mimgneu);
%rates = rates - mean(rates(:));
rates = rates / std(rates(:));
% rates = single(max(rates+2, 0)) + 0.5;

rates = reshape(filter(kernel, 1, rates(:, :)), size(rates));
rates = permute(rates, [2 3 1]);

rates = single(rates);

icell = 0;

xs = repmat(1:Npix, Npix, 1);
ys = xs';

map = zeros(Npix, 'single');
nCoef = .5 + rand(Ncells,1)/2;
nCoef(:) = .8;

meanR = mean(R0(:));
while 1
    % propose new center 
    xcenter(1+icell) = ceil(rand * (Npix-2*radius) + radius);
    ycenter(1+icell) = ceil(rand * (Npix-2*radius) + radius);
    
    ds = (xcenter - xcenter(1+icell)).^2 + (ycenter - ycenter(1+icell)).^2;
    ds = ds(1:icell).^.5;
    
    if icell==0 || min(ds)> radius * 1.33 %2/3
        icell = icell + 1;
        
        dall = (xs - xcenter(icell)).^2 + (ys - ycenter(icell)).^2;
        dall = dall.^.5;
        
        ipix{icell} = find(dall<radius-1e-3);
        npix = numel(ipix);
        
        lam{icell} = 1./(1 + exp(-(dall(ipix{icell}) - 2)));
        lam{icell} = lam{icell}/mean(lam{icell}(:));
        
        iy = ceil(ycenter(icell) *10/Npix);
        ix = ceil(xcenter(icell) *10/Npix);
        
        R =  2*rand * R0(:, iy, ix)/meanR;
%         R = rand * mean(R0(:)) * ones(size(R0,1), 1); 
%         R = mean(R0(:)) * ones(size(R0,1), 1); 
                
        sp = muf * poissrnd(R/10);
        
        spall(:,icell) = sp;
        
        amp  = exprnd(1) * amp0 * mimgneu(1,iy,ix)/mean(mimgneu(:));
        Fcell = amp * filter(kernel, 1, sp(:));
        Fcell = Fcell + .35 * mean(Fcell(:)) * exprnd(1); 
                
        map(ipix{icell}) =  map(ipix{icell}) + 1;
        
        F(:, icell) = Fcell;
        
        if icell==Ncells
            break;
        end
    end
end

F = F(ntpad+1:end, :);
rates = rates(:,:,ntpad+1:end);
R0 = R0(ntpad+1:end, :, :);
spall = spall(ntpad+1:end, :);

mALL = zeros(Npix - pixpad, Npix - pixpad, size(rates,3), 'int16');

nbatch = 1000;
meanImg = zeros(Npix, Npix);

for i = 1:ceil(NT/nbatch)
    trange = (i-1)*nbatch + [1:nbatch];
    trange(trange>NT) = [];
    
    mov = single(imresize(rates(:,:,trange), [Npix Npix]));
    mov = reshape(mov, Npix^2, []);
    
    for icell = 1:Ncells
        mov(ipix{icell}, :) = bsxfun(@plus, nCoef(icell) * mov(ipix{icell}, :), ...
            lam{icell}(:) * F(trange,icell)');
    end
    
    mov = mov .* max(0, 1+ sig_noise*randn(size(mov), 'single'));
    mov = reshape(mov, Npix, Npix, nbatch);
    
    meanImg = meanImg + sum(mov,3);
%     imagesc(meanImg);
%     drawnow;
    
    %movShift = zeros(Npix/5,Npix/5,NT,'single');
%     disp('shifting now');
    
    % shift movie
    for j = 1:nbatch
        for d = 1:2
            mshift      = round(xyshift(trange(j),d));
            indmap{d}  = [max(1,mshift+1) : min(Npix,Npix+mshift)];
            indorig{d} = [max(1,-1*mshift+1) : min(Npix,Npix-mshift)];
        end
        mj = zeros(Npix,Npix,'single');
        mc = mov(:,:,j); % current frame
        mj(indmap{1},indmap{2}) = mc(indorig{1},indorig{2});
        mov(:,:,j)             = mj(1:end,1:end);
    end
    
    mALL(:,:,trange) = 100*mov(pixpad/2+[1:Npix-pixpad],pixpad/2+[1:Npix-pixpad],:);
    
    fprintf('Frame %d shifted \n', i*nbatch);
end

meanImg     = meanImg/NT;
meanImg     = meanImg(pixpad/2+[1:Npix-pixpad],pixpad/2+[1:Npix-pixpad],:);

    
