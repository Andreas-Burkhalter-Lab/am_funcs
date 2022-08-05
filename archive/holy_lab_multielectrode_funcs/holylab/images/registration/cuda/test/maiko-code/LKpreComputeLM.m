function [VTdWdP H Hlm] = LKpreComputeLM(image1,b,xn,jac_x,jac_y,VTdWdP,options)
% precompuation of hessian and steepest descent images
% equation references from:  
% Baker, S., & Matthews, I. (2004). Lucas-Kanade 20 Years On: A Unifying Framework. International Journal of Computer Vision, 56(3), 221-255. doi:10.1023/B:VISI.0000011205.11775.fd

[Gy Gx] = gradient(image1,options.grad_sig,options.grad_sig); % padded gradient
Gx = Gx(b+1:end-b,b+1:end-b); % remove padding
Gy = Gy(b+1:end-b,b+1:end-b); % remove padding
% calculate deepest descent matrix
for j = 1:xn % VT*dW/dp in equations 35 and 36
    VTdWdP(j,:) = [Gx(j) Gy(j)] * [jac_x(j,:) ; jac_y(j,:)];
end
% calculate hessian 
% equation 36
H = VTdWdP(:,:)'*VTdWdP(:,:);
% calculate diagonal hessian for Levenberg-Marquadt
% equation 92 (second half) 
Hlm = zeros(size(H,1));
for i = 1:size(H,1) % calculate diagonal hessian
    Hlm(i,i) = H(i,i);
end

end