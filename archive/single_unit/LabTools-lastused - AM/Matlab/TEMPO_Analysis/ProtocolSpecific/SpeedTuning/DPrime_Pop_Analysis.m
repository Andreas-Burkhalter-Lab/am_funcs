
%load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\dprime_population_Gamma.mat');
%load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\dprime_population_LogGaussFree.mat');
%load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\dprime_population_LogGaussFixed0_3.mat');
load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\dprime_population_LogGaussFixed0_01.mat');
size(d_prime_pop)

%load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\Gammaderiv_population.mat');
%load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\LogGaussFreederiv_population.mat');
%load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\LogGaussFixed0_3deriv_population.mat');
load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\LogGaussFixed0_01deriv_population.mat');
size(deriv_pop)

spd = 0.05 : 0.05: 32.0;

if (length(spd) ~= size(d_prime_pop, 2))
    display('dimensions do not match');
    return;
end

d_prime_mean = nanmean(abs(d_prime_pop));
d_prime_median = nanmedian(abs(d_prime_pop));

d_prime_maunsell = sqrt(nansum(abs(d_prime_pop).^2));  %same as sqrt Fisher information

figure;
%plot(log(spd), log(d_prime_mean), 'r-', log(spd), log(d_prime_median), 'm-');
%plot(log(spd), log(d_prime_mean), 'r-');
hold on;
%plot(log(spd), log(d_prime_maunsell), 'b-');
plot((spd), (d_prime_maunsell), 'b-');
xlabel('speed');
ylabel('d');

figure;
hold on;
%plot(log(spd), log(1./(d_prime_avg)), 'r-');
%plot(log(spd), log(1./(d_prime_maunsell)), 'b-');
%plot((spd), (1./(d_prime_mean)), 'r-', (spd), (1./(d_prime_median)), 'm-');
%plot((spd), (1./(d_prime_mean)), 'r-');

thresh = (1./(d_prime_maunsell));
plot((spd), thresh, 'b-');

xlabel('speed (deg/s)');
ylabel('d-prime=1 threshold (deg/s)');

Weber = thresh./spd;
figure;
plot((spd), Weber, 'b-');

xlabel('speed (deg/s)');
ylabel('Weber fraction');

% 
% raw_slopes = load('Z:\LabTools\Matlab\TEMPO_Analysis\ProtocolSpecific\SpeedTuning\SpeedTuningRawSlopes.dat');
% 
% raw_diffs = raw_slopes;
% %multiply through by delta speed
% raw_diffs(:,1) = raw_slopes(:,1)*0.5;
% raw_diffs(:,2) = raw_slopes(:,2)*0.5;
% raw_diffs(:,3) = raw_slopes(:,3)*1.0;
% raw_diffs(:,4) = raw_slopes(:,4)*2.0;
% raw_diffs(:,5) = raw_slopes(:,5)*4.0;
% raw_diffs(:,6) = raw_slopes(:,6)*8.0;
% raw_diffs(:,7) = raw_slopes(:,7)*16.0;
% 
% [vals, indxs] = max(raw_slopes');
% figure; 
% counts = hist(indxs,6);
% 
% %[vals, indxs] = max(raw_diffs');
% %figure; hist(indxs,6);
