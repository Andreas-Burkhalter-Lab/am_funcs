function[ f2 ] = warpFrame( f, A, T )
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