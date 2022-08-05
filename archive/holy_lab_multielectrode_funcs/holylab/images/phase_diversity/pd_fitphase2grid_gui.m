function A = pd_fitphase2grid_gui(phi,pupildata,gridfunc)
% pd_fitphase2grid_gui: fit the actuator grid to phase-diversity results
%
% Syntax:
%   A = pd_fitphase2grid_gui(phi,pupildata)
%   A = pd_fitphase2grid_gui(phi,pupildata,gridfunc)
% where
%   phi is a structure array, one per actuator. Each element has fields
%     actuator: the number of the given actuator
%     phiV OR Zindex, Zvalue: either the raw phase or the Zernike-fit
%       phase parameters (see PD_ANALYZE_VDEP_STACK)
% and pupildata has fields
%     H0, rho, theta: pupil data
% The optional input gridfunc has default value @mirao52_layout, and must
% have a syntax consistent with that function.
%
% This GUI will ask the user to "flip" the phases so that they have the
% desired sign (phases under an odd-transformation have no impact on the
% resulting image).  Then the user hits return, and clicks on the peaks
% of three actuators.  This then generates an affine transformation that
% maps the actuator positions onto the phase plane.
%
% See also: PD_ANALYZE_VDEP_STACK.
  
% Copyright 2009 by Timothy E. Holy
  
  if (nargin < 3)
    gridfunc = @mirao52_layout;
  end
  
  n_act = length(phi);
  actuators = [phi.actuator];
  % Calculate the raw phases, if need be
  if ~isfield(phi,'phiV')
    uZindex = unique([phi.Zindex]);
    zernv = zernike_values(pupildata.rho,pupildata.theta,uZindex);
    for actIndex = 1:n_act
      phiV = zeros(size(pupildata.rho));
      Zlist = findainb(phi(actIndex).Zindex,uZindex);
      for Zindex = 1:length(Zlist)
        phiV = phiV + phi(actIndex).Zvalue(Zindex) * zernv(:,:,Zlist(Zindex));
      end
      phi(actIndex).phiV = phiV;
    end
  end
  xr = [-1 1] * max(pupildata.rho(:));
  
  % Plot the phases (each actuator in its own axis)
  dims = CalcSubplotDims(n_act);
  hfig_act = figure;
  hax = zeros(1,n_act);
  hctrl = hax;
  for k = 1:n_act
    hax(k) = subplot(dims(1),dims(2),k);
    imagesc(xr,xr,fftshift(phi(k).phiV))
    title(num2str(phi(k).actuator))
    colorbar
    set(hax(k),'Units','pixels');
    pos = get(hax(k),'Position');
    pos(2) = pos(2)+pos(4); pos(4) = 20; pos(3) = 50;
    hctrl(k) = uicontrol('Parent',hfig_act,'Style','checkbox','Position',pos,'String','Flip','Callback',@(src,event) pdfp2g_flip(src,event,hax(k)));
  end
  set([hax hctrl],'Units','normalized')
  
  input('Flip phases as needed and hit return to continue')
  
  % Get the grid of actuators and generate the geometric
  % transform A that puts the actuators at user-defined points
  [act_grid,act_xy] = gridfunc();
  [clickx,clicky,button,hax_click] = ginput_with_axis(3);
  click_ax_index = findainb(hax_click,hax);
  click_actuators = actuators(click_ax_index);
  if (length(unique(click_actuators)) < length(click_actuators))
    error('Must click 3 different axes')
  end
  X = act_xy(click_actuators,:);
  T = cp2tform(X,[clickx clicky],'affine');
  
  % Plot each actuator in its panel
  act_phase = tformfwd(T,act_xy);
  for k = 1:n_act
    axes(hax(k))
    line(act_phase(actuators(k),1),act_phase(actuators(k),2),'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','w')
  end
  
  % Plot the entire grid of actuators
  figure
  pupilimg = fftshift(pupildata.H0);
  pupilimg = repmat(pupilimg,[1 1 3]);
  pupilimg(pupilimg == 0) = 0.8;
  image(xr,xr,pupilimg)
  for k = 1:size(act_phase,1)
    text(act_phase(k,1),act_phase(k,2),num2str(k),'HorizontalAlignment','center','VerticalAlignment','middle');
  end
  phaserange = [min(act_phase); max(act_phase)];
  set(gca,'XLim',phaserange(:,1),'YLim',phaserange(:,2),'Visible','off')
  axis equal
  
  A = T.tdata.T;

%   xymin = min(act_xy);
%   xymax = max(act_xy);
%   mag = diff(xr) / max(xymax-xymin);
%   xyoffset = (xymin+xymax)/2 * mag;
%   A = [mag 0; 0 mag; -xyoffset];
%   
%   hfig_grid = figure;
  
end

function pdfp2g_flip(src,event,hax)
  him = findobj(hax,'type','image');
  imd = get(him,'CData');
  imd = fliplr(imd);
  imd = flipud(imd);
  imd = - imd;
  set(him,'CData',imd);
end
  
  