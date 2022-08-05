function [] = imflow_display3D(varargin)
% imflow_display3D: Function wrapper of imflow_display to visualize uc in three 
% different directions. The current figure will be automically closed when
% done.
%
% usage:
%       imflow_display3D(uc, options.sz_spatial, direc)
% where
%       uc is the warpping displacement vectors,
%       options.sz_spatial is the pixel size of the images,
%       direc can be 1, 2 or 3 which represents x, y or z axis
% For example,
%       imflow_display3D(uc, options.sz_spatial, 3)
% shows the uc displacement on the x-y plane.

% Copyright 2012 Jian Wang and Time Holy

  
  if ~eq(nargin, 3)
    error('The input arguments must be 3');
  end
    
  uc = varargin{1};
  ucSize = varargin{2};  
  funcHandle = register_gui_utilities;
  unorm = funcHandle.upc2mg(uc,ucSize);
  unormSize = size(unorm);
  
  direc = varargin{3};
  if direc < 1 || direc > 3
    error('The 3rd input value must be 1, 2 or 3');
  end
  
  figure;
  i = 1;
  while ~isnan(i)
    clear ucplane;
    ucplane= [];
    switch direc
      case 3
        ucplane(:,:,1) = unorm(:,:,i,1);
        ucplane(:,:,2) = unorm(:,:,i,2);        
      case 2
        ucplane(:,:,1) = unorm(:,i,:,1);
        ucplane(:,:,2) = unorm(:,i,:,3);
      case 1
        ucplane(:,:,1) = unorm(i,:,:,2);
        ucplane(:,:,2) = unorm(i,:,:,3);
    end    
    clf, imflow_display(ucplane);
    title(['Plot of UC: ',num2str(i),'/',num2str(unormSize(direc))]);
    i = keystepper(1:unormSize(direc),i);
  end
  close(gcf);
  
  
  