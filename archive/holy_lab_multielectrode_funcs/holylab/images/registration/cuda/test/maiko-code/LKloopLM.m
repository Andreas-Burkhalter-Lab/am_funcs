function [Pc Fc i_count_c] = LKloopLM(bIms,cIm,VTdWdP,H,Hlm,process_points,Pc,Fc,Pcomp,i_count_c,cpboxcurr,xg,yg,options,cframe)
% affine transformation:
% W(x,p) is:
% [1+P(1)  P(3)    P(5) ;
%  P(2)    1+P(4)  P(6) ;
%  0       0       1   ]
% inverting dP is same as composing the above matrix with dP and taking the
% inverse of the composed matrix, it's written here explicity for speed
%
% composing P and dP is the same as multiplying the matrices together
% however here I explicitly write the formulas for composing P & dP speed.
% equation references from : 
% Baker, S., & Matthews, I. (2004). Lucas-Kanade 20 Years On: A Unifying Framework. International Journal of Computer Vision, 56(3), 221-255. doi:10.1023/B:VISI.0000011205.11775.fd

for i = process_points
    if options.LevenMarq % if we want to use Levenburg-Marquadt Method
        delta = 0.001; % initialize delta
    else 
        delta = 0;
    end
    itercount = 0;
    P = zeros(1,6);
    P(1:4) = Pc(1:4,i); % give previous update parameters for robustness
    dPer = 1; % arbitrary to enter the loop
    % slice the images
    template = bIms(:,:,i);
    % original image is the template (remove padding
    testing = cIm(cpboxcurr(i,2):cpboxcurr(i,2)+cpboxcurr(i,4),cpboxcurr(i,1):cpboxcurr(i,1)+cpboxcurr(i,3));    
    % precompute the first warp
    % create transformation matrix -- equation 2
    Hg = [(1+P(1)) P(3) P(5) ; P(2) 1+P(4) P(6) ; 0 0 1]; % our guess of the homography (transformation)
    % warp the current image
    warped = affine_transform_2d_double(testing,xg,yg,Hg,false);    
    % calculate the error between pixels
    ierror = warped-template;
    ImEr = sum(ierror(:).^2);
    
    while dPer > options.error_threshold && itercount <  options.max_iter
        itercount = itercount+1; % increment the counter
        % solve for dP (non-inverted)
        % H = (H(:,:,i) + delta * Hlm(:,:,i))
        % use backslash operator rather than inverse of H for speed
        % equation 93
        dP = (H(:,:,i) + delta * Hlm(:,:,i))\(VTdWdP(:,:,i)'*ierror(:));
        % invert dP (equations 32 & 33)
        dP = ((1+dP(1))*(1+dP(4))-dP(2)*dP(3))^(-1)*[-dP(1)-dP(1)*dP(4)+dP(2)*dP(3)  -dP(2)  -dP(3)  -dP(4)-dP(1)*dP(4)+dP(2)*dP(3)  -dP(5)-dP(4)*dP(5)+dP(3)*dP(6)  -dP(6)-dP(1)*dP(6)+dP(2)*dP(5)];
        % enforce some parameters for robustness
        if itercount < options.nonTranslationIterations
            %dP = [0     0    0    0     dP(5) dP(6)];
            dP = [dP(1) 0    0    dP(4) 0     0];
        end
        % compose inverted dP with P
        P_lm = compose_P(P,dP); % equation 16
        % create transformation matrix
        Hg = [(1+P_lm(1)) P_lm(3) P_lm(5) ; P_lm(2) 1+P_lm(4) P_lm(6) ; 0 0 1]; % our guess of the homography (transformation)
        % warp the current image
        warped = affine_transform_2d_double(testing,xg,yg,Hg,false);
        % calculate the error between pixels
        ierror_lm = warped-template;
        ImEr_lm = sum(ierror_lm(:).^2);
        if ImEr < ImEr_lm && options.LevenMarq
            % bad update 
            delta = delta * 10;
        else
            % good update (or not using Levenberg-Marquardt) 
            delta = delta / 10;
            P = P_lm; % update the parameters
            dPer = sum(abs(dP)); % calculate the new epsilon in dP
            ImEr = ImEr_lm; % use the new error
            ierror = ierror_lm; % new error image
        end
        
        % if we want to show the deformation of each region live 
        if options.plotLive && cframe == 3 && ~options.use_par_cpu && options.verbose;
            subplot(1,3,1); imagesc(template)
            title(['itercount : ' num2str(itercount)])
            subplot(1,3,2); imagesc(warped);
            title(['dP epsilon : ' num2str(dPer)])
            subplot(1,3,3); imagesc(ierror); colormap gray
            title(['err : ' num2str(ImEr)])
            drawnow;
            pause(.1)
        end
    end
    i_count_c(i) = itercount;
    Pc(:,i) =  P; % insert the P into our P current 
    Pf =  compose_P(P,Pcomp(:,i)); % compose P with any previous updates for calculating F
    % build the matrix then invert it to find F 
    % we invert it becuase F is relative to the template, whereas the
    % inverse compositional method finds the deformation relative to the
    % image
    res = inv([(1+Pf(1)) Pf(3) Pf(5) ; Pf(2) 1+Pf(4) Pf(6) ; 0 0 1]')'; % build an invert the transform to compute F
    Fc(:,i) = res(:); % insert the result into our F matrix (vectorized for convenience)
end

end


function P = compose_P(P,dP) 
% This composition is the same as building a warp for P and dP:
% W(x,P):                           W(x,dP):
% [1+P(1)  P(3)    P(5) ;           [1+dP(1)  dP(3)    dP(5) ;
%  P(2)    1+P(4)  P(6) ;            dP(2)    1+dP(4)  dP(6) ;
%  0       0       1   ];            0       0       1      ];
% from the parameters P and dP, then muliplying the resulting matrices:
% W(x,P)*W(x,dP)
%
% derived from equation 16
% 
% Here it is calculated explicity for speed: 
P = [P(1)+dP(1)+P(1)*dP(1)+P(3)*dP(2) ...
    P(2)+dP(2)+P(2)*dP(1)+P(4)*dP(2) ...
    P(3)+dP(3)+P(1)*dP(3)+P(3)*dP(4) ...
    P(4)+dP(4)+P(2)*dP(3)+P(4)*dP(4) ...
    P(5)+dP(5)+P(1)*dP(5)+P(3)*dP(6) ...
    P(6)+dP(6)+P(2)*dP(5)+P(4)*dP(6)  ];
end
