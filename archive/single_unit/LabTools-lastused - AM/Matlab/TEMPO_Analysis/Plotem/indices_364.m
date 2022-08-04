% indices364.m defines all the I_xxxx and SIN_xxxx indices for Anal364

% indices to columns of the spike index (*.index) file (index.dat):
 SIN_FIX = 1; SIN_TRGON = 2;
 SIN_STON = 3; SIN_STOFF = 4;
 SIN_SAC = 5; SIN_TRGAC = 6;
 SIN_MSTIM_ON = 7;
 SIN_MSTIM_OFF = 8;
 SIN_TOT = 9;
 SIN_FPON = 10 ; 	 % FPON is always zero and will be the last column
                         % of the index file.  The columns of zeros is
			 % added in loadfiles.m
			 
% indices that are not used in this paradigm
 SIN_TRGOFF = [];
 SIN_FPOFF = []; 
 SIN_REW = [];
 
% indices to column of spikes *.s file (foo.dat).
I_CTR_COH = 1; 
I_CTR_DIR = 2; 
I_CTR_DIAM = 3; 
I_CTR_SPEED = 4; 
I_CTR_HOR_DISP = 5; 
I_CTR_VER_DISP = 6; 
I_CTR_BIN_CORR = 7; 
I_SURR_COH = 8; 
I_SURR_DIR = 9; 
I_SURR_DIAM = 10; 
I_SURR_SPEED = 11; 
I_SURR_HOR_DISP = 12; 
I_SURR_VER_DISP = 13; 
I_SURR_BIN_CORR = 14;
I_SHOW_SURR = 15;
I_MSTIM = 16;
I_FP2FIX = 17; 
I_FP2FIX_T = 18;
I_FIX2ST = 19; 
I_FIX2ST_T = 20;
I_ST = 21; 
I_ST_T = 22; 
I_ST2SAC = 23; 
I_ST2SAC_T = 24;
I_STOFF2SAC = 25; 
I_STOFF2SAC_T = 26;
I_SAC = 27; 
I_SAC_T = 28;
I_MSTIM_PER = 29;
I_MSTIM_PER_T = 30;
I_TRGCHOICE = 31; 
I_CORRECT = 32; 

% indices that are not used in this paradigm
 I_SID = []; I_T1X = []; I_T1Y = [];
 I_T2X = []; I_T2Y = []; I_DIAM = [];


% These are the events that are defined for this paradigm
% we need to declare their names and the mapping between
% them and the columns of the index file.  GDLH 1/10/95
tmp = [...
 'FixPt On|','FixPt Acquire|',...
 'Targ On|','Stimulus On|', 'Stimulus Off|',...
 'Saccade|', 'Targ Acquire'];
trigPopupString = tmp(:)';
trigIndexMap = [...
    SIN_FPON SIN_FIX SIN_TRGON SIN_STON SIN_STOFF SIN_SAC SIN_TRGAC];
eventIndexMap = [nan trigIndexMap];

% protocol codes
DIR_TUNING =			1;		% std. direction tuning run 
HORIZ_DISP_TUNING =	2;		% std. horiz disparity tuning run (w/ monoc controls)
VERT_DISP_TUNING =	3;		% vertical disparity tuning run (w/ monoc controls)
AREA_SUMMATION =		4;		% std. area summation test
SPEED_TUNING =			5;		% std. speed tuning run
DIREC_HDISP_TUNING =	6;		% 2D Expt: Direction vs. H. Disparity (w/ monoc controls)
HDISP_WWO_SURROUND =	7;		% Horiz. disp tuning with/without surround dots
CTR_VS_SURR_HDISP =	8;		% 2D Expt: Center HDisp vs. Surround HDisp
HDISP_VS_VDISP	= 		9;		% 2D Expt: HDisp vs. VDisp (of Center) 
STD_DIREC_DISCRIM =	10;	% Standard Direction Discrim., at multiple disparities, if wanted 
DIREC_DISCRIM_MSTIM =11;	% Same as above, but with microstim. conditions 
DEPTH_DISCRIM_MSTIM =12;	% Depth Discrimination Expt. w/ microstim.
MULTI_HDISP_TUNING  =13;	% Disparity tuning: moving, 0 coh, and static dots
