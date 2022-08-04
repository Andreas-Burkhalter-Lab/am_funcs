%before running, load the PSTHxxx.mat file into the workspace
Path_Defs

MU_range = [1:50];
SU_range = [51:63];
N_smooth_pts=20;

MU_Mov = nansum(mov_PSTH_out(MU_range, :)-spont_PSTH_out(MU_range, :));
MU_Stat = nansum(stat_PSTH_out(MU_range, :)-spont_PSTH_out(MU_range, :));
%MU_Spont = nansum(spont_PSTH_out(MU_range, :));

SU_Mov = nansum(mov_PSTH_out(SU_range, :)-spont_PSTH_out(SU_range, :));
SU_Stat = nansum(stat_PSTH_out(SU_range, :)-spont_PSTH_out(SU_range, :));
%SU_Spont = nansum(spont_PSTH_out(SU_range, :));

MU_Mov2 = BoxcarFilter(MU_Mov,N_smooth_pts);
MU_Stat2 = BoxcarFilter(MU_Stat,N_smooth_pts);
Time2 = BoxcarFilter(time_out,N_smooth_pts);
SU_Mov2 = BoxcarFilter(SU_Mov,N_smooth_pts);
SU_Stat2 = BoxcarFilter(SU_Stat,N_smooth_pts);

%normalize
max_MU = max(MU_Mov2);
max_SU = max(SU_Mov2);
MU_Mov2 = MU_Mov2/max_MU;
MU_Stat2 = MU_Stat2/max_MU;
SU_Mov2 = SU_Mov2/max_SU;
SU_Stat2 = SU_Stat2/max_SU;

figure;
plot(Time2-30, [MU_Mov2' MU_Stat2']);
figure;
plot(Time2-30, [SU_Mov2' SU_Stat2']);

Out = [Time2' MU_Mov2' MU_Stat2' SU_Mov2' SU_Stat2'];
outfile = 'Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\HDispTuning\Static_Moving_PSTH_out.dat';
save(outfile, 'Out', '-ASCII');