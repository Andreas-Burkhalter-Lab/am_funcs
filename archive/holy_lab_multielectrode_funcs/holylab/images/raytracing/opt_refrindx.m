function n = opt_refrindx(glass,wavelength)
% n = opt_refrindx(glass,wavelength)
% 
% OPT_REFRINDX refractive index calculations.
% GLASS - name of glass type (a string), WAVELENGTH in nanometers.
% n: the index of refraction of that glass at that wavelength
%
% Alternatively you can make 'glass' a number, and that number is returned
% as n.
%
% Data from JML Optical Industries, Inc available at:
%   http://www.jmloptical.com/level2/index.aspx?pagename=TechInfo/transmissiongraphs.aspx
% and from Thorlabs: (for c0550, tac4, and fds9)
%   http://www.thorlabs.com/Equation.cfm?Section=2&Ref=12
% and from
%   http://www.luxpop.com/
% for water & saline
%   http://www.corning.com/docs/specialtymaterials/pisheets/H0607_CaF2_Prod
%   uct_Sheet.pdf
% for CaF2

% Copyright B. Gustavsson 20050804

persistent glass_names_sellmeier refr_consts_sellmeier

if isempty(glass_names_sellmeier)
  
  qwe = '';
  fp = fopen('Sellmeier.glass.refr','r');
  while ~feof(fp)
    qwe = str2mat(qwe,fgetl(fp));
  end
  glass_names_sellmeier = lower(qwe(4:end,1:7));
  refr_consts_sellmeier = str2num(qwe(4:end,8:end));
end

lambda_ref = [5876 4861 6563]*1e-1; % Input wavelength in nm

I = strmatch(lower(glass),glass_names_sellmeier,'exact');
if ~isempty(I)
  lambda = wavelength*1e-3; % change units to micrometer
  A = refr_consts_sellmeier(I,1:3);
  B = refr_consts_sellmeier(I,4:6);
  n = sqrt(1+A(1)*lambda.^2./(lambda.^2-B(1)) + A(2)*lambda.^2./(lambda.^2-B(2))+A(3)*lambda.^2./(lambda.^2-B(3)));
  return
end

switch lower(glass)
 case 'air'
  nref = [1 1 1];
 case 'water'
  %nref = [1.33386 1.33797 1.33224];
  %nref = [1.33387 1.33797 1.33198];
  lref = [238 245 275 313 365 407 435 546 589 633 780 900 1300 1550];
  nref = [1.38541 1.38107 1.36736 1.35668 1.34792 1.3434 1.34113 1.3353 ...
	  1.33383 1.33257 1.32943 1.32744 1.32109 1.31585];
  n = interp1(lref,nref,wavelength);
  return
 case 'saline'
  lref = [407 435 546 589 633];
  nref = [1.346 1.343 1.337 1.336 1.335];      % For 10g salt/L, approx 170mM total salt
  n = interp1(lref,nref,wavelength);
  return
  %nref = [1.336 1.340 1.334];   % For 10g salt/L, approx 170mM total salt
 case 'cytop'
  lref = [238 245 275 313 365 407 435 546 589 633 1300 1550];
  nref = [1.35764 1.35637 1.35393 1.35132 1.34840 1.34566 1.34404 1.3402 ...
	  1.34 1.3395 1.3348 1.3335];
  n = interp1(lref,nref,wavelength);
  return 
 case 'caf2'
  lref = [250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000];
  nref = [1.4673 1.45399 1.44652 1.44185 1.43872 1.43649 1.43485 1.43358 1.43259 1.43178 1.43112 1.43055 1.43006 1.42963 1.42925 1.4289];
  n = interp1(lref,nref,wavelength);
  return
 case 'acrylic'
  nref = [1.491 1.496 1.488];
 case 'b270'
  nref = [1.5230 1.5292 1.5202];
 case 'baf10'
   nref = [1.67002 1.67999 1.66579];
 case 'bak1'
  nref = [1.5725 1.5794 1.5694];
 case 'bak2'
  nref = [1.5399 1.5462 1.5372];
 case 'bak4'
  nref = [1.56883 1.5759 1.56576];
 case 'balkn3'
  nref = [1.51849 1.52447 1.51586];
 case 'bk7'
  nref = [1.5168 1.5224 1.5143];
 case 'f2'
  nref = [1.6200 1.6320 1.6150];
 case 'f3'
  nref = [1.61293 1.62461 1.60806];
 case 'f4'
  nref = [1.6165 1.6284 1.6116];
 case 'fusedsilica'
  nref = [1.458 1.463 1.456];
 case 'k5'
  nref = [1.5224 1.5285 1.5198];
 case 'k7'
  nref = [1.51112 1.517 1.50854];
 case 'lasfn9'
  nref = [1.850 1.8689 1.8425];
 case 'lah71'
  nref = [1.8502 1.8689 1.8425];
 case 'pyrex'
  nref = [1.473 1.478 1.471];
 case 'sapphire'
  nref = [1.7682 1.7756 1.7649];
 case 'sf1'
  nref = [1.71736 1.73462 1.71031];
 case 'sf2'
  nref = [1.6476 1.6612 1.6421];
 case 'sf5'
  nref = [1.6727 1.6875 1.66661];
 case 'sf8'
  nref = [1.6889 1.7046 1.6825];
 case 'sf10'
  nref = [1.72825 1.74648 1.72085];
 case 'sf11'
  nref = [1.7847 1.8064 1.7759];
 case 'sf12'
  nref = [1.64831 1.66187 1.64271];
 case 'sf15'
  nref = [1.69895 1.71546 1.69221];
 case 'sf18'
  nref = [1.7215 1.7390 1.7143];
 case 'sf19'
  nref = [1.6668 1.6811 1.6609];
 case 'sf56'
  nref = [1.7847 1.8061 1.7760];
 case 'sk3'
  nref = [1.6088 1.6160 1.6056];
 case 'sk5'
  nref = [1.5891 1.5958 1.5861];
 case 'sk11'
  nref = [1.5638 1.5702 1.5610];
 case 'sk16'
  nref = [1.6204 1.6275 1.6172];
 case 'ssk2'
  nref = [1.6223 1.63048 1.61878];
 case 'ssk4a'
  nref = [1.61765 1.62547 1.61427];
 case 'ssk51'
  nref = [1.60361 1.61147 1.60022];
 case 'zk5'
  nref = [1.53375 1.54049 1.53084];
 case 'c0550'
  nref = [1.60500 1.61341 1.60141];
 case 'tak4'
  nref = [1.73400 1.74403 1.72965];
 case 'fds9'
  nref = [1.84666 1.87204 1.83651];
 case 'h-k9l'
  C1 = 2.2702566; C2 = -9.19881011e-3; C3 = 1.16097061e-2; C4 = -7.61239111e-5; C5 = 2.85587271e-5; C6 = -1.25664861e-6;
  lamda = wavelength * 1e-3;
  n = sqrt(C1 + C2 * lamda.^2 + C3 * lamda.^-2  + C4 * lamda.^-4 + C5 * lamda.^-6 + C6 * lamda.^-8);
  return
 case 'h-lak2'
  C1 = 2.8081523; C2 = -1.4159869e-2; C3 = 1.7693707e-2; C4 = 7.6788759e-4; C5 = -7.4886174e-5; C6 = 4.7044756e-6;
  lamda = wavelength * 1e-3;
  n = sqrt(C1 + C2 * lamda.^2 + C3 * lamda.^-2  + C4 * lamda.^-4 + C5 * lamda.^-6 + C6 * lamda.^-8);
  return 
 case 'h-zf6'
  C1 = 2.96731720E+00; C2 = -1.3989564E-02; C3 = 3.34566800E-02; C4 = 3.0717804E-03; C5 = -2.6682434E-04; C6 = 2.89502410E-05;
  lamda = wavelength * 1e-3;
  n = sqrt(C1 + C2 * lamda.^2 + C3 * lamda.^-2  + C4 * lamda.^-4 + C5 * lamda.^-6 + C6 * lamda.^-8);
  return
 case 'h-zf52'
  C1 = 3.2606321; C2 = -1.7780398E-02; C3 = 4.0902938E-02; C4 = 5.6076934E-03; C5 = -5.6434039E-04; C6 = 5.4763391E-5;
  lamda = wavelength * 1e-3;
  n = sqrt(C1 + C2 * lamda.^2 + C3 * lamda.^-2  + C4 * lamda.^-4 + C5 * lamda.^-6 + C6 * lamda.^-8);
  return  
 case 'h-lak59'
  C1 = 2.7821195; C2 = 3.8698879e-3; C3 = 3.2299384e-2; C4 = -3.1791589e-3; C5 = 4.3583685e-4; C6 = -2.06224134e-5;
  lamda = wavelength * 1e-3;
  n = sqrt(C1 + C2 * lamda.^2 + C3 * lamda.^-2  + C4 * lamda.^-4 + C5 * lamda.^-6 + C6 * lamda.^-8);
  return  
 case 'h-zf7l'
  C1 = 3.117148300E+000; C2 = -1.161627700E-002; C3 = 4.486497500E-002; C4 = 1.658849100E-003; C5 = 2.571685400E-005; C6 = 1.475765700E-005
  lamda = wavelength * 1e-3;
  n = sqrt(C1 + C2 * lamda.^2 + C3 * lamda.^-2  + C4 * lamda.^-4 + C5 * lamda.^-6 + C6 * lamda.^-8);
  return  
 case 'zf7l'
  C1 = 3.121597700E+000; C2 = -1.187173600E-002; C3 = 4.008587300E-002; C4 = 3.498813900E-003; C5 = -2.584040300E-004; C6 = 2.841529000E-005;
  lamda = wavelength * 1e-3;
  n = sqrt(C1 + C2 * lamda.^2 + C3 * lamda.^-2  + C4 * lamda.^-4 + C5 * lamda.^-6 + C6 * lamda.^-8);
  return  
 case 'n138'
  n = 1.38;
  return;
 case 'n140'
  n = 1.4;
  return;
 % From here on out, we only have index at the d-line (587.6nm) and the
 % V-number (the ratio (n_d - 1) / (n_f - n_c), with the f-line at 486.1nm
 % and the c-line at 656.3nm.
 % The first number is n_d, the second the V-number
 % So, these will be pretty inaccurate outside of the visible
  case 'fk51'
    nref = [1.48656 84.47];
  case 'pcd4'
    nref = [1.61800 63.39];
  case 'sf3'
    nref = [1.74000 28.24];
  case 'laf28'
    nref = [1.77250 46.96];  % This one seems to vary widely!
  case 'k3'
    nref = [1.518 59.0];    
 otherwise
  n = str2double(glass); % Try a number encoded as a string
  if isempty(n) || isnan(n)
    error(['Material "' glass '" not recognized']);
  else
    return
  end
  %nref = [1 1 1];
end


% If we only have n_d and the V-number, just calculate it linearly:
if (length(nref) == 2)
  n_d = nref(1);
  V = nref(2);
  dlambda_ref = -diff(lambda_ref(2:3));
  dlambda = wavelength - lambda_ref(1);
  m = (n_d - 1)/n_d/V/dlambda_ref;
  n = n_d * (1 + m*dlambda);
  return
end
%normx = (1./lambda_ref.^2 - 1/wavelength^2)*wavelength^2;

%[abc,s,mu] = polyfit(1./lambda_ref.^2,nref,2);
%n = polyval(abc,1/wavelength^2,s,mu);

xr = 1./lambda_ref.^2;
xlambda = 1./wavelength.^2;
% Use Neville's method to evaluate polynomial
for j = 1:length(xlambda)
  x = xlambda(j);
  for i = 1:2
    p(i) = ((x-xr(i+1))*nref(i) + (xr(i)-x)*nref(i+1))./(xr(i)-xr(i+1));
  end
  n(j) = ((x-xr(3))*p(1) + (xr(1)-x)*p(2))./(xr(1)-xr(3));
end
