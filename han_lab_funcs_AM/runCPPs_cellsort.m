close all; clear all;
saveDIR='G:\MA Data\Interneurons\';
saveName='CPPstats.mat';
if ~exist(saveDIR,'dir')
    mkdir(saveDIR)
end
    varlist='mouseCPP_CS';
m=1;
curr=1;
[mouseCPP_CS(m).FPyrs{curr},mouseCPP_CS(m).FInts{curr},mouseCPP_CS(m).Rotations{curr},mouseCPP_CS(m).Forwards{curr},...
    mouseCPP_CS(m).CorrsPyr{curr},mouseCPP_CS(m).CorrsInts{curr},mouseCPP_CS(m).DistsPyr{curr},mouseCPP_CS(m).DistsInt{curr},mouseCPP_CS(m).CPPPyr{curr},mouseCPP_CS(m).CPPInt{curr},mouseCPP_CS(m).ybinned{curr},...
    mouseCPP_CS(m).rewards{curr},...
    mouseCPP_CS(m).CorrsRunPyr{curr},mouseCPP_CS(m).CorrsRunInt{curr},mouseCPP_CS(m).CorrswoezPyr{curr},mouseCPP_CS(m).CorrswoezInt{curr},...
    mouseCPP_CS(m).mfballout{curr},mouseCPP_CS(m).mtballout{curr},mouseCPP_CS(m).mtaballout{curr},mouseCPP_CS(m).mvrout{curr},mouseCPP_CS(m).mvraout{curr},mouseCPP_CS(m).mPRTAoutPyr{curr},mouseCPP_CS(m).mPRTAoutInt{curr},...
    mouseCPP_CS(m).dFcorrPyr{curr},mouseCPP_CS(m).dFcorrInt{curr},mouseCPP_CS(m).PyrActivity{curr},mouseCPP_CS(m).IntActivity{curr}]=...
    MixedNeuronTestCPP(5.2,3,'S21 Pre Cellsort',...
    {'170330_SW_000_002__XCNC_1_cellsort_F.mat','170330_SW_000_002__XCNC_2_cellsort_F.mat',...
    '170330_SW_000_002__XCNC_3_ca_F.mat'},...
    {'G:\CPP\S21\Pre\','G:\CPP\S21\Pre\',...
    'G:\CPP\S21\Pre\'},0,'new');
save([saveDIR,saveName],varlist);

m=1;
curr=2;
[mouseCPP_CS(m).FPyrs{curr},mouseCPP_CS(m).FInts{curr},mouseCPP_CS(m).Rotations{curr},mouseCPP_CS(m).Forwards{curr},...
    mouseCPP_CS(m).CorrsPyr{curr},mouseCPP_CS(m).CorrsInts{curr},mouseCPP_CS(m).DistsPyr{curr},mouseCPP_CS(m).DistsInt{curr},mouseCPP_CS(m).CPPPyr{curr},mouseCPP_CS(m).CPPInt{curr},mouseCPP_CS(m).ybinned{curr},...
    mouseCPP_CS(m).rewards{curr},...
    mouseCPP_CS(m).CorrsRunPyr{curr},mouseCPP_CS(m).CorrsRunInt{curr},mouseCPP_CS(m).CorrswoezPyr{curr},mouseCPP_CS(m).CorrswoezInt{curr},...
    mouseCPP_CS(m).mfballout{curr},mouseCPP_CS(m).mtballout{curr},mouseCPP_CS(m).mtaballout{curr},mouseCPP_CS(m).mvrout{curr},mouseCPP_CS(m).mvraout{curr},mouseCPP_CS(m).mPRTAoutPyr{curr},mouseCPP_CS(m).mPRTAoutInt{curr},...
    mouseCPP_CS(m).dFcorrPyr{curr},mouseCPP_CS(m).dFcorrInt{curr},mouseCPP_CS(m).PyrActivity{curr},mouseCPP_CS(m).IntActivity{curr}]=...
    MixedNeuronTestCPP(5.2,3,'S21 Post Cellsort',...
    {'170417_SW_000_001_NCaccnew1_cellsort_F.mat','170417_SW_000_001__XC_2_cellsort_F.mat',...
    '170417_SW_000_001__XC_3_cellsort_F.mat'},...
    {'G:\CPP\S21\Post\','G:\CPP\S21\Post\',...
    'G:\CPP\S21\Post\'},0,'new');
save([saveDIR,saveName],varlist);

m=2;
curr=1;
[mouseCPP_CS(m).FPyrs{curr},mouseCPP_CS(m).FInts{curr},mouseCPP_CS(m).Rotations{curr},mouseCPP_CS(m).Forwards{curr},...
    mouseCPP_CS(m).CorrsPyr{curr},mouseCPP_CS(m).CorrsInts{curr},mouseCPP_CS(m).DistsPyr{curr},mouseCPP_CS(m).DistsInt{curr},mouseCPP_CS(m).CPPPyr{curr},mouseCPP_CS(m).CPPInt{curr},mouseCPP_CS(m).ybinned{curr},...
    mouseCPP_CS(m).rewards{curr},...
    mouseCPP_CS(m).CorrsRunPyr{curr},mouseCPP_CS(m).CorrsRunInt{curr},mouseCPP_CS(m).CorrswoezPyr{curr},mouseCPP_CS(m).CorrswoezInt{curr},...
    mouseCPP_CS(m).mfballout{curr},mouseCPP_CS(m).mtballout{curr},mouseCPP_CS(m).mtaballout{curr},mouseCPP_CS(m).mvrout{curr},mouseCPP_CS(m).mvraout{curr},mouseCPP_CS(m).mPRTAoutPyr{curr},mouseCPP_CS(m).mPRTAoutInt{curr},...
    mouseCPP_CS(m).dFcorrPyr{curr},mouseCPP_CS(m).dFcorrInt{curr},mouseCPP_CS(m).PyrActivity{curr},mouseCPP_CS(m).IntActivity{curr}]=...
     MixedNeuronTestCPP(5.2,4,'S22 Pre Cellsort',...
    {'170330_SW_000_003__XCNC_1_cellsort_F.mat','170330_SW_000_003__XCNC_2_cellsort_F.mat',...
    '170330_SW_000_003__XCNC_3_cellsort_F.mat','170330_SW_000_003__XCNC_4_cellsort_F.mat',...
    },...
    {'G:\CPP\S22\Pre\','G:\CPP\S22\Pre\',...
    'G:\CPP\S22\Pre\','G:\CPP\S22\Pre\',...
    },0,'new');
save([saveDIR,saveName],varlist);

m=2;
curr=2;
[mouseCPP_CS(m).FPyrs{curr},mouseCPP_CS(m).FInts{curr},mouseCPP_CS(m).Rotations{curr},mouseCPP_CS(m).Forwards{curr},...
    mouseCPP_CS(m).CorrsPyr{curr},mouseCPP_CS(m).CorrsInts{curr},mouseCPP_CS(m).DistsPyr{curr},mouseCPP_CS(m).DistsInt{curr},mouseCPP_CS(m).CPPPyr{curr},mouseCPP_CS(m).CPPInt{curr},mouseCPP_CS(m).ybinned{curr},...
    mouseCPP_CS(m).rewards{curr},...
    mouseCPP_CS(m).CorrsRunPyr{curr},mouseCPP_CS(m).CorrsRunInt{curr},mouseCPP_CS(m).CorrswoezPyr{curr},mouseCPP_CS(m).CorrswoezInt{curr},...
    mouseCPP_CS(m).mfballout{curr},mouseCPP_CS(m).mtballout{curr},mouseCPP_CS(m).mtaballout{curr},mouseCPP_CS(m).mvrout{curr},mouseCPP_CS(m).mvraout{curr},mouseCPP_CS(m).mPRTAoutPyr{curr},mouseCPP_CS(m).mPRTAoutInt{curr},...
    mouseCPP_CS(m).dFcorrPyr{curr},mouseCPP_CS(m).dFcorrInt{curr},mouseCPP_CS(m).PyrActivity{curr},mouseCPP_CS(m).IntActivity{curr}]=...
    MixedNeuronTestCPP(5.2,5,'S22 Post Cellsort',...
    {'170417_SW_000_002__XCNC_1_cellsort_F.mat','170417_SW_000_002__XC_2_cellsort_F.mat',...
    '170417_SW_000_002__XC_3_cellsort_F.mat','170417_SW_000_002__XC_4_cellsort_F.mat',...
    '170417_SW_000_002__XC_5_cellsort_F.mat'},...
    {'G:\CPP\S22\Post\','G:\CPP\S22\Post\',...
    'G:\CPP\S22\Post\','G:\CPP\S22\Post\',...
    'G:\CPP\S22\Post\'},0,'new');
save([saveDIR,saveName],varlist);
