function CPP_plots(mouseCPP,m,curr,name)
pname='G:\MA Data\CPPs\Figures\';
if ~exist(pname,'dir')
    mkdir(pname);
end
if ~exist([pname,'\',name,'\'],'dir')
    mkdir([pname,'\',name,'\']);
end
import mlreportgen.dom.*;
figs=Document(strrep([pname,name,'Report'],'.',''),'html-file');
% figs=Document(strrep([pname,name,'Report'],'.',''),'pdf');

allStuff=Table;

%%
maxh=876.2818;
useinds=1:round(5.2*60*15);
yb=(((mouseCPP(m).ybinned{curr}{1}(useinds,1))/maxh)+1)*1.5;
down=yb<1;
up=yb>2;
% if ~isempty(mouseCPP(m).spks{curr})
% Sp=mouseCPP(m).spks{curr}{1}(useinds,:);
% end
% Fc2=mouseCPP(m).FPyrs{curr}{1};
F=calc_Fc3(mouseCPP(m).FPyrs{curr}{1}(useinds,:));
forward=mouseCPP(m).Forwards{curr}{1}(useinds,:);
rotation=mouseCPP(m).Rotations{curr}{1}(useinds,:);
%% Cell Rois on Mean Image Stacks
maskPics=TableRow;
disp_rois_S2P_html(mouseCPP(m).meanImages{curr},mouseCPP(m).masks{curr},maskPics,pname,name);
append(allStuff,maskPics);

%% Specificity Info
specFigs=TableRow;
plot_Prefs_cpp(yb,up,down,F,pname,name,'Fc3',specFigs);
plot_speed_preference(F,forward,rotation,pname,name,'Ball Speed',specFigs);
plot_speed_preference_vr(F,[0; diff(yb(:,1))],pname,name,'VR Speed',specFigs);
append(allStuff,specFigs);



%% Look at place activity lets keep ends in for now
%with mid
midPlaces=TableRow;
[fc,MIs,prior]=find_places_CPP(F,yb,36);
plot_places(tuning_curve_processing(fc),[],prior,pname,name,'all',midPlaces);
%w/o mid
use=yb<1 | yb>2;
[fc,MIs,prior]=find_places_CPP(F(use,:),yb(use),36);
plot_places(tuning_curve_processing(fc),[],prior,pname,name,'noMid',midPlaces);

angleRow=TableRow;
plot_angle(F,up,down,pname,name,angleRow)
append(allStuff,midPlaces);
append(allStuff,angleRow);

%% ensemble angle

% angle=dot(mFdo,mFup)/(norm(mFdo)*norm(mFup));

%% Malahanobis distance mean pop vector 1s bins


close all;
append(figs,allStuff);
close(figs)

end