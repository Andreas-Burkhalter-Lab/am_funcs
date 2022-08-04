function depth = depth_from_disp(disp, fix_depth, inter_oc)
%-----------------------------------------------------------------------
%function depth = depth_from_disp(disp, fix_depth, inter_oc)
%calculate depth from disparity
%inputs are:
%disp = disparity of object
%fix_depth = viewing distance
%inter_oc = inter ocular distance
%formula is:
%(d((N*d)/a))/(1-((N*d)/a))
%-----------------------------------------------------------------------
%side1 = atan(inter_oc/(2*fix_depth)) - disp/2
%side2 = atan(inter_oc/(2*(fix_depth+0)))
%depth = 3;
disp = disp *(pi/180);
num = inter_oc;
den = 2*(tan(atan(inter_oc/(2*fix_depth)) - (disp/2)));

delta_d = (num/den) - fix_depth;
depth = fix_depth + delta_d;