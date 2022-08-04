% The function videostabilize takes as input: (1) a data structure frames with field im that contains the original video sequence;
%                                             (2) the number of Gaussian pyramid levels L;
%                                              3) a binary image roi that specifies the region of interest, with value of 1 for those pixels to be considered in the motion estimation, and 0 otherwise. 
% It is assumed that each frame is a grayscale image. The output of this function is a data structure stable with field im that contains the stabilized video sequence, and a data structure motion
% with fields A and T that contain the estimated affine and translation parameters.
%http://www.cs.dartmouth.edu/farid/downloads/publications/tr07.pdf
% -------------------------------------------------------------------------
%%% STABILIZE VIDEO
function[ Tout, Aout, stable ] = videostabilize( frames, roi, L )
N = size(frames,3);
roiorig = roi;
%% ESTIMATE PAIRWISE MOTION
Acum = [1 0 ; 0 1];
Tcum = [0 ; 0];
% stable(1).roi = roiorig;
affine(N-1,:,:)=zeros(2);
translation(N-1,:)=[0;0];
for k = 1 : N-1
% [A,T] = opticalflow(frames(:,:,k+1),frames(:,:,k),roi,L);%( frames(k+1).im, frames(k).im, roi, L ); 
[A,T] = opticalflow(frames(:,:,k+1),frames(:,:,k),roi,L,k);
affine(k,:,:) = A;
translation(k,:) = T;
[Acum,Tcum] = accumulatewarp( Acum, Tcum, A, T );
roi = warp( roiorig, Acum, Tcum );
Tout(:,k)=Tcum;
Aout(:,:,k)=Acum;
end
%% STABILIZE TO LAST FRAME
 stable(:,:,N)=frames(:,:,N);%stable(N).im =frames(:,:,N); % frames(N).im;
Acum = [1 0 ; 0 1];
Tcum = [0 ; 0];
for k = N-1 : -1 : 1
[Acum,Tcum] = accumulatewarp( Acum, Tcum, squeeze(affine(k,:,:)), translation(k,:)');
stable(:,:,k)=warp(frames(:,:,k),Acum,Tcum); %stable(k).im = warp(frames(k).im, Acum, Tcum );
end
% -------------------------------------------------------------------------
%% ALIGN TWO FRAMES (f2 to f1)
function[ Acum, Tcum ] = opticalflow( f1, f2, roi, L, frame )
f2orig = f2;
Acum = [1 0 ; 0 1];
Tcum = [0 ; 0];
for k = L : -1 : 0
%%% DOWN-SAMPLE
f1d = down( f1, k );
f2d = down( f2, k );
ROI = down( roi, k );
%%% COMPUTE MOTION
[Fx,Fy,Ft] = spacetimederiv( f1d, f2d );
[A,T] = computemotion( Fx, Fy, Ft, ROI, frame, k );
T = (2^k)*T;
%A=A2A1 T=A2T1+T2
[Acum,Tcum] = accumulatewarp( Acum, Tcum, A, T );
%%% WARP ACCORDING TO ESTIMATED MOTION
f2 = warp( f2orig, Acum, Tcum );
% f2 = warp( f2, A, T);
end
% -------------------------------------------------------------------------
%% COMPUTE MOTION
function[ A, T ] = computemotion( fx, fy, ft, roi, frame, level )
[ydim,xdim] = size(fx);
[x,y] = meshgrid( [1:xdim]-xdim/2, [1:ydim]-ydim/2 );
%% TRIM EDGES
fx   = fx( 3:end-2, 3:end-2 );
fy   = fy( 3:end-2, 3:end-2 );
ft   = ft( 3:end-2, 3:end-2 );
roi  = roi( 3:end-2, 3:end-2 );
x    = x( 3:end-2, 3:end-2 );
y    = y( 3:end-2, 3:end-2 );
ind  = find( roi > 0 );
x    = x(ind); y = y(ind);
fx   = fx(ind); fy = fy(ind); ft = ft(ind);
xfx  = x.*fx; xfy = x.*fy; yfx  = y.*fx; yfy = y.*fy;
% A = [m1,m2;m3,m4] T=[2^(l-1)m5;2^(l-1)m6]
% Error minimize E=[ft - (m1x+m2y+m5-x)fx -(m3x+m4y+m6-y)fy]^2 = [k-cTm]^2
% where k= ft+xfx+yfy, cT=(xfx yfx xfy yfy fx fy)
%Derivate and set to zero, get: m=inv(ccT)(ck)
M(1,1) = sum( xfx .*xfx );  M(1,2) = sum( xfx .*yfx );  M(1,3) = sum( xfx .*xfy );
M(1,4) = sum( xfx .*yfy );  M(1,5) = sum( xfx .*fx );   M(1,6) = sum( xfx .*fy );
M(2,1) = M(1,2);            M(2,2) = sum( yfx .*yfx );  M(2,3) = sum( yfx .* xfy );
M(2,4) = sum( yfx .*yfy );  M(2,5) = sum( yfx .*fx );   M(2,6) = sum( yfx .*fy );
M(3,1) = M(1,3);            M(3,2) = M(2,3);            M(3,3) = sum( xfy .*xfy );
M(3,4) = sum( xfy .*yfy );  M(3,5) = sum( xfy .*fx );   M(3,6) = sum( xfy .*fy );
M(4,1) = M(1,4);            M(4,2) = M(2,4);            M(4,3) = M(3,4);
M(4,4) = sum( yfy .*yfy );  M(4,5) = sum( yfy .*fx ); 	M(4,6) = sum( yfy .*fy );
M(5,1) = M(1,5);            M(5,2) = M(2,5);            M(5,3) = M(3,5);
M(5,4) = M(4,5);            M(5,5) = sum( fx .*fx );    M(5,6) = sum( fx .*fy );
M(6,1) = M(1,6);            M(6,2) = M(2,6);            M(6,3) = M(3,6);
M(6,4) = M(4,6);            M(6,5) = M(5,6);            M(6,6) = sum( fy .*fy );

k = ft + xfx + yfy;

b(1) = sum( k .*xfx );     b(2) = sum( k .*yfx );
b(3) = sum( k .*xfy );     b(4) = sum( k .*yfy );
b(5) = sum( k .*fx );      b(6) = sum( k .*fy );

% sing_test=rcond(M);

% warning('error','MATLAB:singularMatrix');
% 
% try
%     v = M\b'; %mwa
%     % v = inv(M)*b';
% catch
%     disp(['Found a singular, frame ', num2str(frame), ' level ', num2str(level)]);
%     M
%     rcond(M)
%     v = pinv(M)*b' %mwa 2
% end

if rcond(M)>(1e-6)
    v=M\b';
elseif sum(M(:))==0
    v=[1 0 0 1 0 0];
    frame
else
    v=pinv(M)*b';
    frame
end
A = [v(1) v(2) ; v(3) v(4)];
T = [v(5) ; v(6)];

% -------------------------------------------------------------------------
%% WARP IMAGE
function[ f2 ] = warp( f, A, T )
[ydim,xdim] = size( f );
[xramp,yramp] = meshgrid( (1:xdim)-xdim/2, (1:ydim)-ydim/2 );
P       = [xramp(:)' ; yramp(:)'];
% A=[1,0;0,1];
% %Apply affine transformation
P       = A*P;
xramp2  = reshape( P(1,:), ydim, xdim ) + T(1);
yramp2  = reshape( P(2,:), ydim, xdim ) + T(2);
% xramp2=xramp+T(1);
% yramp2=yramp+T(2);
f2      = interp2( xramp, yramp, f, xramp2, yramp2, 'cubic' ); % warp 
ind     = find( isnan(f2) );
f2(ind) = 0;
% -------------------------------------------------------------------------
%%% BLUR AND DOWNSAMPLE (L times)
function[ f ]  = down( f, L )
blur = [1 2 1]/4;
for k = 1 : L
f = conv2( conv2( f, blur, 'same' ), blur', 'same' );
f = f(1:2:end,1:2:end);
end
% -------------------------------------------------------------------------
%%% SPACE/TIME DERIVATIVES
function[ fx, fy, ft ] = spacetimederiv( f1, f2 )
%%% DERIVATIVE FILTERS
pre   = [0.5 0.5];
deriv = [0.5 -0.5];
%%% SPACE/TIME DERIVATIVES
fpt = pre(1)*f1 + pre(2)*f2; % pre-filter in time; avg f1 and f2
fdt = deriv(1)*f1 + deriv(2)*f2; % differentiate in time; difference b/t f1 and f2
fx = conv2( conv2( fpt, pre', 'same' ), deriv, 'same' );  
fy = conv2( conv2( fpt, pre, 'same' ), deriv', 'same' );
ft = conv2( conv2( fdt, pre', 'same' ), pre, 'same' );
% -------------------------------------------------------------------------
%%% ACCUMULATE WARPS
function[ A2, T2 ] = accumulatewarp( Acum, Tcum, A, T )
A2 = A*Acum;
T2 = A*Tcum + T;



