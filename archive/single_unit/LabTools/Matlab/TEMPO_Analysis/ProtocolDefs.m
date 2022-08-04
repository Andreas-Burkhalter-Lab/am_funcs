%----------------------------------------------------------------------------------------------------------
%-- ProtocolDefs.m: Contains defines for protocol-specific analyses.  
%--	GCD, 1/3/2000
%-- Last revised 3/1/2002 BJP
%----------------------------------------------------------------------------------------------------------

%defines for various pre-canned experimental protocols
%NOTE: these defines MUST match those in the TEMPO protocol, (ProtDefs.PRO)
NUM_DOTS_PROTOCOLS = 22;
DIRECTION_TUNING = 0; 
SPEED_TUNING = 1;
HDISP_TUNING = 2;
DEPTH_DISCRIM = 3;
DIRECTION_DISCRIM = 4;
HDISP_GRAD_TUNING = 5;
SIZE_TUNING = 6;
SURROUND_MAPPING = 7;
SURROUND_TUNING = 8;
RELATIVE_DISPARITY = 9;
STEREOACUITY = 10;
EYE_CALIBRATION = 11;
RF_MAPPING = 12;
TRANS_RELATIVE_DISPARITY = 13;
SIM_DIST_DISP_ONLY = 14;
DEPTH_DISCRIM_NOVAR = 15;
SIM_DIST_VERG_ONLY = 16;
SIM_DIST_DISP_VERG = 17;
SIM_DIST_CURVATURE_DISCRIM = 18;
CUED_DIRECTION_DISCRIM = 19;
AXIS_CUED_DIRECTION_TUNING = 20;
ORIENT_CUE_DIRECTION_TUNING = 21;

%----------------------------------------------------
% Grating protocols
%----------------------------------------------------
NUM_GRATING_PROTOCOLS = 7;
GRATING_ORIENTATION = 22;
GRATING_SPATIAL_FREQ = 23;
GRATING_TEMPORAL_FREQ = 24;
GRATING_CONTRAST = 25;
GRATING_SIZE = 26;
MOTION_COHERENCE = 27;  %Note: had to make this name different from the protocol name in .log file due to naming conflict with MOOG parameter
GRATING_RF_MAP = 28;

%Defines for MOOG protocols, added by GCD 6/27/03
NUM_MOOG_PROTOCOLS = 16;
DIRECTION_TUNING_3D = 100;
HEADING_DISCRIM = 101;
DIR3D_VESTIB_ACCEL = 103;
DIR3D_VARY_FIXATION = 104;
DIR2D_CUE_CONFLICT = 105;
AZIMUTH_TUNING_1D = 107;
AZ_1D_VARY_FIXATION = 108;
EL_1D_VARY_FIXATION = 109;
HEADING_DISCRIM_FIXONLY = 110;
ROTATION_TUNING_3D = 112;
TILT_TRANSLATION=113; %added by AHC
RVOR_PURSUIT = 114;

% Defines for REVCORR protocols.
NUM_REVCORR_PROTOCOLS = 1;
DIRECTION_REVCORR = 106;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BINDING PROTOCOLS
%% ADDED BY BJP 1/3/01
BIND_DIR_DISC = 200;
BIND_BACKCOLOR = 201;
BIND_BACKHDISP = 202;
BIND_OBJ_LUM = 203;
BIND_OBJ_WIDTH = 204;
BIND_AP_RLUM = 205;
BIND_OBJ_POSITION = 206;
BIND_DOT_DENSITY = 207;
BIND_COHER_DISC = 208;
BIND_CLOSURE_DISC = 209;

FIXATION = 250;
FIX_1_23_45 = 251;
FIX_1_23 = 252;
BIND_EYE_CALIBRATION = 253;
BIND_DIR_TUNING = 254;
BIND_SPATFREQ_TUNING = 255;
BIND_TEMPFREQ_TUNING = 256;
BIND_RF_MAPPING = 257;
BIND_HDISP_TUNING = 258;
BIND_BAR_DIR_TUNING = 259;
BIND_BAR_SPEED_TUNING = 260;
BIND_BAR_HDISP_TUNING = 261;
FIX_VARY_HISTORY = 262;
FIX_VARY_BACKCOLOR = 263;
FIX_TRANS_QUAD = 264;

NUM_BINDING_PROTOCOLS = 24;

%--------------------------------------------------
%Surface orientation protocols - JDN 01/22/04
%--------------------------------------------------
NUM_SURF_PROTOCOLS = 7;
SURF_DEPTH_TUNING = 300;
SURF_SPEED_TUNING = 301;
SURF_RF_MAPPING = 302;
SURF_DIRECTION_TUNING = 303;
SURF_SIZE_TUNING = 304;
SURF_TUNING = 305;
SURF_EYE_CALIBRATION = 306;
SLANT_TUNING = 307;

%--------------------------------------------------
%Motion parallax protocols - JWN 09/16/04
%--------------------------------------------------
NUM_MP_PROTOCOLS = 1;
MOTION_PARALLAX_FIX = 400;

%Protocol description strings
%NOTE: the order of these strings must match the order of the indices above, and there can be NO GAPS
%in the range of indices, GCD
protocol_names = 	{'Direction Tuning',
    'Speed Tuning',
    'Horizontal Disparity Tuning',
    'Depth Discrimination',
    'Direction Discrimination',
    'Horiz. Disparity Gradient Tuning',
    'Size Tuning',
    'Surround Mapping',
    'Surround Tuning',
    'Relative Disparity',
    'Stereoacuity',
    'Eye Calibration',
    'RF Mapping',
    'Transparent Relative Disparity',
    'Simulate Distance Disparities Only',
    'Depth Discrim. Var/Novar',
    'Simulate Distance Vergence Only',
    'Simulate Distance Disparity+Vergence',
    'Simulate Distance Curvature Discrim',
    'Cued Direction Discrimination',
    'Axis Cued Direction Tuning',
    'Orientation Cued Direction Tuning',
    };

%added by GCD 6/27/03 for MOOG protocols
protocol_names(101:106) = {
    'MOOG: 3D Direction Tuning',
    'MOOG: Heading Discrimination', 
    'Undefined',
    'MOOG: Direc 3D Vestib. vary Accel',
    'MOOG: 3D Direction Tuning Gaze'
    'MOOG: 2D Direction Tuning conflict'
};

protocol_names(107) = {
    'REVCORR: Reverse Correlation'
};

protocol_names(108) = {
    'MOOG: AZIMUTH_TUNING_1D'
};

protocol_names(109) = {
    'MOOG: 1D Az Direction Tuning Gaze'
};

protocol_names(110) = {
    'MOOG: 1D El Direction Tuning Gaze'
};
protocol_names(111) = {
    'MOOG: Heading Fixation'
};
protocol_names(112) = {
    'Underfined'
};
protocol_names(113) = {
    'MOOG: ROTATION_TUNNG_3D'
};
 protocol_names(115) = {
    'MOOG: RVOR_PURSUIT'   
};
protocol_names(114) = {
    'MOOG: TILT_TRANSLATION'   
}; % added by AHC

%%% added by BJP 1/3/01
protocol_names(201:210) = {'Binding: Std Dir Disc',
    'Binding: Vary Background Color',
    'Binding: Vary Background Hor. Disp',
    'Binding: Vary Obj. Lum',
    'Binding: Vary Obj. Width',
    'Binding: Vary Aperture Red Lum.',
    'Binding: Vary Obj. Position'
    'Binding: Vary Bkgnd Dot Density'
    'Binding: Coherence Disc' 
    'Binding: Closure Disc'     
};
protocol_names(251:265) = {'Binding: Fixation',
    'Binding: Fix Obj 1 / 23/ 45',
    'Binding: Fix Obj 1 / 23',
    'Eye Calibration',
    'Direction Tuning',
    'Spatial Freq Tuning',
    'Temporal Freq Tuning',
    'RF Mapping',
    'H Disp Tuning',
    'Direction Tuning (Bars)',
    'Speed Tuning (Bars)',
    'H Disp Tuning (Bars)',
    'Binding: Vary History',
    'Binding: Vary Bkgnd Color',
    'Fix Trans Quad'
};

protocol_names(301:308) = {'Depth Tuning Curve',
'Speed Tuning Curve',
'RF Map', 
'Direction Tuning Curve',
'Size Tuning Curve', 
'Surface Orientation Tuning Curves',
'Eye Calibration',
'Tilt by Slant Tuning Curves'};

%Motion parallax protocols - JWN 09/16/04
protocol_names(MOTION_PARALLAX_FIX + 1) = {'Motion Parallax Fix'};

% Grating protocols
protocol_names(GRATING_ORIENTATION+1:GRATING_ORIENTATION+NUM_GRATING_PROTOCOLS) = ...
    {'Grating Orientation',
     'Grating Spatial Freq',
     'Grating Temporal Freq',
     'Grating Contrast',
     'Grating Size',
     'Motion Coherence',  %not actually a grating protocol but put here for numbering convenience, GCD
     'Grating RF Map',
     };

%defines for patch indices
PATCH1 = 1;
PATCH2 = 2;
PATCH3 = 3;
PATCH4 = 4;
NUM_PATCHES = 4;

%defines for MOOG movement parameters
MOOG = 1;
CAMERAS = 2;
NUM_MOOG_ITEMS = 2;  %check this against .log file

%defines for target indices; these are indices into data.targ_params
FP = 1;
T1 = 2;
T2 = 3;
T3 = 4;
T4 = 5;
NUM_TARGETS = 5;

%other defines for binding protocol
NUM_OBJECTS = 6;
NUM_BARS = 12;

%define for cue protocol
NUM_CUES = 1;

%stimtype defines for surface orientation protocol - JDN - 1/22/04
CONGRUENT = 0;
SPEED_ONLY = 1;
DISPARITY_ONLY = 2;
TEXTURE_ONLY = 3;

%constants that indicate control conditions in h. disparity tuning runs
LEYE_CONTROL = -99.0;
REYE_CONTROL = 99.0;
UNCORR_CONTROL = 98.0;

%constant that define cue type in cued paradigms - VR 9/22/05
CUEONLY = 2;
VALID = 1;
NEUTRAL = 0;
INVALID = -1;
NOCUE = -2;

%here, I set up some lists of keywords that MATLAB will understand when reading
%in the .log file from TEMPO, which will contain all parameter values for
%each trial in an experiment, GCD 12/29/99
%First, a list of keywords for dots parameters
num_keys = 1;
dots_keywords{num_keys}='DOTS_DIREC';				DOTS_DIREC=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_SPEED';				DOTS_SPEED=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_COHER';				DOTS_COHER=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_HDISP';				DOTS_HDISP=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_VDISP';				DOTS_VDISP=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_BIN_CORR';			DOTS_BIN_CORR=num_keys; 			num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_HGRAD_MAG';			DOTS_HGRAD_MAG=num_keys; 			num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_HGRAD_ANGLE';		DOTS_HGRAD_ANGLE=num_keys; 		num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_COH_TYPE';			DOTS_COH_TYPE=num_keys; 			num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_DENSITY';			DOTS_DENSITY=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_DWELL';				DOTS_DWELL=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_RLUM';				DOTS_RLUM=num_keys; 					num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_GLUM';				DOTS_GLUM=num_keys; 					num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_BLUM';				DOTS_BLUM=num_keys; 					num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_CONTRAST';			DOTS_CONTRAST=num_keys; 			num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_XCTR';			DOTS_AP_XCTR=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_YCTR';			DOTS_AP_YCTR=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_XSIZ';			DOTS_AP_XSIZ=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_YSIZ';			DOTS_AP_YSIZ=num_keys; 				num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_OFF_ANG';		DOTS_AP_OFF_ANG=num_keys; 			num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_OFF_RAD';		DOTS_AP_OFF_RAD=num_keys; 			num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_STATUS';			DOTS_AP_STATUS=num_keys; 			num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_SHAPE';			DOTS_AP_SHAPE=num_keys; 			num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_AP_DEPTH_ORDER';	DOTS_AP_DEPTH_ORDER=num_keys; 	num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_BIN_CORR_SEED';	DOTS_BIN_CORR_SEED=num_keys; 		num_keys=num_keys+1;
dots_keywords{num_keys}='NOVAR_FLAG';	    NOVAR_FLAG=num_keys; 		num_keys=num_keys+1;
% Additional parameters from depth simulation MLM 7/2/2002
dots_keywords{num_keys}='DEPTH_SETTING';        DEPTH_SETTING=num_keys;             num_keys=num_keys+1;
dots_keywords{num_keys}='DEPTH_DIST_SIM';       DEPTH_DIST_SIM=num_keys;            num_keys=num_keys+1;
dots_keywords{num_keys}='DEPTH_FIX_SIM';        DEPTH_FIX_SIM=num_keys;             num_keys=num_keys+1;
dots_keywords{num_keys}='DEPTH_FIX_REAL';       DEPTH_FIX_REAL=num_keys;            num_keys=num_keys+1;
dots_keywords{num_keys}='DEPTH_SIM_CONFLICT';   DEPTH_SIM_CONFLICT=num_keys;        num_keys=num_keys+1;

%surface orientation params
dots_keywords{num_keys}='DOTS_STIM_TYPE';       DOTS_STIM_TYPE=num_keys;          num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_XSIZE';           DOTS_XSIZE=num_keys;              num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_YSIZE';           DOTS_YSIZE=num_keys;              num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_XSIZE_DISP';      DOTS_XSIZE_DISP=num_keys;              num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_YSIZE_DISP';      DOTS_YSIZE_DISP=num_keys;              num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_POLYSIZE';        DOTS_POLYSIZE=num_keys;           num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_SLANT';           DOTS_SLANT=num_keys;              num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_TILT';            DOTS_TILT=num_keys;               num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_SCRN_LEFT';       DOTS_SCRN_LEFT=num_keys;          num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_SCRN_RIGHT';      DOTS_SCRN_RIGHT=num_keys;         num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_SCRN_TOP';        DOTS_SCRN_TOP=num_keys;           num_keys=num_keys+1;
dots_keywords{num_keys}='DOTS_SCRN_BOTTOM';     DOTS_SCRN_BOTTOM=num_keys;        num_keys=num_keys+1;


NUM_DOTS_PARAMS = num_keys - 1; %use last number


%List of Keywords for MOOG paradigm parameters, added 6/27/03 GCD
num_keys = 1;
moog_keywords{num_keys}='AZIMUTH';				AZIMUTH=num_keys; 				num_keys=num_keys+1;
moog_keywords{num_keys}='ELEVATION';			ELEVATION=num_keys; 			num_keys=num_keys+1;
moog_keywords{num_keys}='AMPLITUDE';			AMPLITUDE=num_keys; 			num_keys=num_keys+1;
moog_keywords{num_keys}='ORIGIN_X';			    ORIGIN_X=num_keys; 			    num_keys=num_keys+1;
moog_keywords{num_keys}='ORIGIN_Y';			    ORIGIN_Y=num_keys; 			    num_keys=num_keys+1;
moog_keywords{num_keys}='ORIGIN_Z';			    ORIGIN_Z=num_keys; 			    num_keys=num_keys+1;
moog_keywords{num_keys}='NUM_SIGMAS';			NUM_SIGMAS=num_keys; 			num_keys=num_keys+1;
moog_keywords{num_keys}='DURATION';			    DURATION=num_keys; 			    num_keys=num_keys+1;
moog_keywords{num_keys}='STIM_TYPE';			STIM_TYPE=num_keys; 			num_keys=num_keys+1;
moog_keywords{num_keys}='HEADING';			    HEADING=num_keys; 			    num_keys=num_keys+1;
moog_keywords{num_keys}='LIFETIME';			    LIFETIME=num_keys; 			    num_keys=num_keys+1;
moog_keywords{num_keys}='LIFETIME_STATUS';		LIFETIME_STATUS=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='COHERENCE';			COHERENCE=num_keys; 			num_keys=num_keys+1;
moog_keywords{num_keys}='FIX_X';			    FIX_X=num_keys; 			    num_keys=num_keys+1;
moog_keywords{num_keys}='FIX_Y';			    FIX_Y=num_keys; 			    num_keys=num_keys+1;
moog_keywords{num_keys}='ROT_DURATION';			ROT_DURATION=num_keys; 			num_keys=num_keys+1;
moog_keywords{num_keys}='ROT_ELEVATION';	    ROT_ELEVATION=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='ROT_AZIMUTH';			ROT_AZIMUTH=num_keys; 			num_keys=num_keys+1;
moog_keywords{num_keys}='ROT_AMPLITUDE';		ROT_AMPLITUDE=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='ROT_NUM_SIGMAS';		ROT_NUM_SIGMAS=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='INTEROCULAR_DIST';		INTEROCULAR_DIST=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_LEYE_RLUM';		STARS_LEYE_RLUM=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_LEYE_GLUM';		STARS_LEYE_GLUM=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_LEYE_BLUM';		STARS_LEYE_BLUM=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_REYE_RLUM';		STARS_REYE_RLUM=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_REYE_GLUM';		STARS_REYE_GLUM=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_REYE_BLUM';		STARS_REYE_BLUM=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_LUM_MULT';		STARS_LUM_MULT=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_DENSITY';		STARS_DENSITY=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_SIZE_X';		    STARS_SIZE_X=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_SIZE_Y';		    STARS_SIZE_Y=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_VOLUME_X';		STARS_VOLUME_X=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_VOLUME_Y';		STARS_VOLUME_Y=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_VOLUME_Z';		STARS_VOLUME_Z=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='STARS_TYPE';		    STARS_TYPE=num_keys; 		    num_keys=num_keys+1;
% More keywords added for Rotation 08/09/05
moog_keywords{num_keys}='GL_CENTER_X';		    GL_CENTER_X=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='GL_CENTER_Y';		    GL_CENTER_Y=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='GL_CENTER_Z';		    GL_CENTER_Z=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='TT_MODE';		        TT_MODE=num_keys; 		        num_keys=num_keys+1;
moog_keywords{num_keys}='FP_ROTATE';		    FP_ROTATE=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='NUM_X_PATCHES';	    NUM_X_PATCHES=num_keys; 	    num_keys=num_keys+1;
moog_keywords{num_keys}='NUM_Y_PATCHES';	    NUM_Y_PATCHES=num_keys; 	    num_keys=num_keys+1;
moog_keywords{num_keys}='PATCH_DWELL';		    PATCH_DWELL=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='PATCH_DENSITY';	    PATCH_DENSITY=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='HEAD_CENTER_X';		HEAD_CENTER_X=num_keys;         num_keys=num_keys+1;
moog_keywords{num_keys}='HEAD_CENTER_Y';		HEAD_CENTER_Y=num_keys; 	    num_keys=num_keys+1;
moog_keywords{num_keys}='HEAD_CENTER_Z';		HEAD_CENTER_Z=num_keys; 	    num_keys=num_keys+1;
% More keywords added by JWN 01/31/05
% Not all are actually used.
moog_keywords{num_keys}='MP_TRIAL_TYPE';		MP_TRIAL_TYPE=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='PATCH_DEPTH';		    PATCH_DEPTH=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='MOVEMENT_PHASE';		MOVEMENT_PHASE=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='MOVE_MAGNITUDE';		MOVE_MAGNITUDE=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='MOVE_EXPONENT';		MOVE_EXPONENT=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='MOVE_SIGMA';		    MOVE_SIGMA=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='MOVE_FREQUENCY';		MOVE_FREQUENCY=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='FIELDSIZE_X';		    FIELDSIZE_X=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='FIELDSIZE_Y';		    FIELDSIZE_Y=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='FIELDLOC_X';		    FIELDLOC_X=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='FIELDLOC_Y';		    FIELDLOC_Y=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='BACKGROUND_DENSITY';	BACKGROUND_DENSITY=num_keys;	num_keys=num_keys+1;
moog_keywords{num_keys}='CLIP_NEAR';		    CLIP_NEAR=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='CLIP_FAR';		        CLIP_FAR=num_keys; 		        num_keys=num_keys+1;
moog_keywords{num_keys}='FIXATION_SIZE';		FIXATION_SIZE=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='PATCH_COHERENCE';		PATCH_COHERENCE=num_keys; 		num_keys=num_keys+1;
moog_keywords{num_keys}='DOTS_DISP';		    DOTS_DISP=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='PATCH_SCALE';		    PATCH_SCALE=num_keys; 		    num_keys=num_keys+1;
moog_keywords{num_keys}='PATCH_DRIFT';		    PATCH_DRIFT=num_keys; 		    num_keys=num_keys+1;
NUM_MOOG_PARAMS = num_keys - 1;

num_keys = 1;
revcorr_keywords{num_keys} = 'PATCH_DIMS';      PATCH_DIMS = num_keys;           num_keys = num_keys + 1;
revcorr_keywords{num_keys} = 'DOT_DWELL';       DOT_DWELL = num_keys;            num_keys = num_keys + 1;
revcorr_keywords{num_keys} = 'DOT_DENSITY';     DOT_DENSITY = num_keys;          num_keys = num_keys + 1;
revcorr_keywords{num_keys} = 'DOT_DIREC';       DOT_DIREC = num_keys;            num_keys = num_keys + 1;
revcorr_keywords{num_keys} = 'DOT_SPEED';       DOT_SPEED = num_keys;            num_keys = num_keys + 1;
revcorr_keywords{num_keys} = 'DOT_DISPARITY';   DOT_DISPARITY = num_keys;        num_keys = num_keys + 1;
revcorr_keywords{num_keys} = 'CORRDELAY_LOW';   CORRDELAY_LOW = num_keys;        num_keys = num_keys + 1;
revcorr_keywords{num_keys} = 'CORRDELAY_HIGH';  CORRDELAY_HIGH = num_keys;       num_keys = num_keys + 1;
revcorr_keywords{num_keys} = 'CORRDELAY_INC';   CORRDELAY_INC = num_keys;        num_keys = num_keys + 1;
NUM_REVCORR_PARAMS = num_keys - 1;

num_keys = 1;
gratings_keywords{num_keys} = 'GRAT_ORIENTATION';   GRAT_ORIENTATION = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRAT_SPATIAL_FREQ';   GRAT_SPATIAL_FREQ = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRAT_TEMPORAL_FREQ';   GRAT_TEMPORAL_FREQ = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'DRAW_GRATINGS';   DRAW_GRATINGS = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRATINGS_INC';   GRATINGS_INC = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRATINGS_AP';   GRATINGS_AP = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRAT_WIDTH';   GRAT_WIDTH = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRAT_HEIGHT';   GRAT_HEIGHT = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRAT_XCTR';   GRAT_XCTR = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRAT_YCTR';   GRAT_YCTR = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRAT_CONTRAST';   GRAT_CONTRAST = num_keys; num_keys = num_keys + 1;
gratings_keywords{num_keys} = 'GRAT_TYPE';   GRAT_TYPE = num_keys; num_keys = num_keys + 1;
NUM_GRATINGS_PARAMS = num_keys - 1;

%Now, a list of keywords for target parameters
num_keys = 1;
targ_keywords{num_keys}='TARG_XCTR';			TARG_XCTR=num_keys; 			num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_YCTR';			TARG_YCTR=num_keys; 			num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_ZCTR';			TARG_ZCTR=num_keys; 			num_keys=num_keys+1;  %added for MOOG protocols, 6/27/03 GCD
targ_keywords{num_keys}='TARG_XSIZ';			TARG_XSIZ=num_keys; 			num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_YSIZ';			TARG_YSIZ=num_keys; 			num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_ZSIZ';			TARG_ZSIZ=num_keys; 			num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_SHAPE';			TARG_SHAPE=num_keys;            num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_HDISP';			TARG_HDISP=num_keys;        	num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_RLUM';			TARG_RLUM=num_keys; 			num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_GLUM';			TARG_GLUM=num_keys; 			num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_BLUM';			TARG_BLUM=num_keys; 			num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_LUM_MULT';		TARG_LUM_MULT=num_keys;     	num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_WIN_FWIDTH';      TARG_WIN_FWIDTH=num_keys;   	num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_WIN_FHEIGHT';     TARG_WIN_FHEIGHT=num_keys;      num_keys=num_keys+1;
targ_keywords{num_keys}='TARG_PATCH_SIZE';      TARG_PATCH_SIZE=num_keys;       num_keys=num_keys+1;
NUM_TARG_PARAMS = num_keys - 1; %use last number

%Now, a list of cue parameters, VR 04/01/05
num_keys = 1;
cue_keywords{num_keys}='CUE_DIREC';             CUE_DIREC=num_keys;             num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_XCTR';              CUE_XCTR=num_keys;              num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_YCTR';              CUE_YCTR=num_keys;              num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_XSIZ';              CUE_XSIZ=num_keys;              num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_YSIZ';              CUE_YSIZ=num_keys;              num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_HDISP';             CUE_HDISP=num_keys;             num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_RLUM';              CUE_RLUM=num_keys;              num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_GLUM';              CUE_GLUM=num_keys;              num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_BLUM';              CUE_BLUM=num_keys;              num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_LUM_MULT';          CUE_LUM_MULT=num_keys;          num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_STATUS';            CUE_STATUS=num_keys;            num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_STIM_ISI';          CUE_STIM_ISI = num_keys;        num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_VALIDITY';          CUE_VALIDITY = num_keys;        num_keys=num_keys+1;
cue_keywords{num_keys}='CUE_TYPE';              CUE_TYPE=num_keys;              num_keys=num_keys+1;
cue_keywords{num_keys}='AXIS_CUE_STATUS';       AXIS_CUE_STATUS = num_keys;     num_keys=num_keys+1;
cue_keywords{num_keys}='AXIS_CUE_DIREC';        AXIS_CUE_DIREC = num_keys;      num_keys=num_keys+1;
cue_keywords{num_keys}='DIREC_SPREAD';          DIREC_SPREAD = num_keys;        num_keys=num_keys+1;
NUM_CUE_PARAMS = num_keys - 1;

%Now, a list of keywords for miscellaneous parameters that need to get dropped *EVERY* trial.  
num_keys = 1;
misc_keywords{num_keys}='OUTCOME';					OUTCOME=num_keys; 				num_keys=num_keys+1;
misc_keywords{num_keys}='MICROSTIM';				MICROSTIM=num_keys; 			num_keys=num_keys+1;
misc_keywords{num_keys}='CONDITION';				CONDITION=num_keys; 			num_keys=num_keys+1;
NUM_MISC_PARAMS = num_keys - 1;

num_keys = 1;
neuron_keywords{num_keys}='PREFERRED_SPATFREQ';			PREFERRED_SPATFREQ=num_keys; 	num_keys=num_keys+1;
%% indices for Pref_dir to RF diam must match those in one_time_params
neuron_keywords{num_keys}='PREFERRED_DIRECTION';	    PREFERRED_DIRECTION=num_keys;   num_keys=num_keys+1;
neuron_keywords{num_keys}='PREFERRED_SPEED';			PREFERRED_SPEED=num_keys; 		num_keys=num_keys+1;
neuron_keywords{num_keys}='PREFERRED_HDISP';			PREFERRED_HDISP=num_keys; 		num_keys=num_keys+1;
neuron_keywords{num_keys}='PREFERRED_ORIENTATION';	    PREFERRED_ORIENTATION=num_keys;   num_keys=num_keys+1;
neuron_keywords{num_keys}='PREFERRED_SPATIAL_FREQ';		PREFERRED_SPATIAL_FREQ=num_keys; 		num_keys=num_keys+1;
neuron_keywords{num_keys}='PREFERRED_TEMPORAL_FREQ';	PREFERRED_TEMPORAL_FREQ=num_keys; 		num_keys=num_keys+1;
neuron_keywords{num_keys}='PREFERRED_DEPTH';			PREFERRED_DEPTH=num_keys; 		num_keys=num_keys+1;
neuron_keywords{num_keys}='RF_XCTR';					RF_XCTR=num_keys; 				num_keys=num_keys+1;
neuron_keywords{num_keys}='RF_YCTR';					RF_YCTR=num_keys; 				num_keys=num_keys+1;
neuron_keywords{num_keys}='RF_DIAMETER';				RF_DIAMETER=num_keys; 			num_keys=num_keys+1;
neuron_keywords{num_keys}='PREFERRED_TEMPFREQ';			PREFERRED_TEMPFREQ=num_keys; 	num_keys=num_keys+1;

%Now, a list of keywords for miscellaneous parameters that only get dropped on the *FIRST* trial.  
num_keys = 1;
one_time_keywords{num_keys}='METHOD';					METHOD=num_keys; 				num_keys=num_keys+1;
one_time_keywords{num_keys}='PREFERRED_DIRECTION';	    PREFERRED_DIRECTION=num_keys;   num_keys=num_keys+1;
one_time_keywords{num_keys}='PREFERRED_SPEED';		    PREFERRED_SPEED=num_keys; 		num_keys=num_keys+1;
one_time_keywords{num_keys}='PREFERRED_AZIMUTH';		PREFERRED_AZIMUTH=num_keys;     num_keys=num_keys+1; %added for MOOG, 6/27/03 GCD
one_time_keywords{num_keys}='PREFERRED_ELEVATION';		PREFERRED_ELEVATION=num_keys;   num_keys=num_keys+1; %added for MOOG, 6/27/03 GCD
one_time_keywords{num_keys}='PREFERRED_HDISP';		    PREFERRED_HDISP=num_keys; 		num_keys=num_keys+1;
one_time_keywords{num_keys}='PREFERRED_ORIENTATION';	PREFERRED_ORIENTATION=num_keys;   num_keys=num_keys+1;
one_time_keywords{num_keys}='PREFERRED_SPATIAL_FREQ';	PREFERRED_SPATIAL_FREQ=num_keys;   num_keys=num_keys+1;
one_time_keywords{num_keys}='PREFERRED_TEMPORAL_FREQ';	PREFERRED_TEMPORAL_FREQ=num_keys;   num_keys=num_keys+1;
one_time_keywords{num_keys}='RF_XCTR';					RF_XCTR=num_keys; 				num_keys=num_keys+1;
one_time_keywords{num_keys}='RF_YCTR';					RF_YCTR=num_keys; 				num_keys=num_keys+1;
one_time_keywords{num_keys}='RF_DIAMETER';			    RF_DIAMETER=num_keys; 			num_keys=num_keys+1;
one_time_keywords{num_keys}='VERG_WIN_FWIDTH';		    VERG_WIN_FWIDTH=num_keys; 		num_keys=num_keys+1;
one_time_keywords{num_keys}='VERG_WIN_FHEIGHT';		    VERG_WIN_FHEIGHT=num_keys; 	    num_keys=num_keys+1;
one_time_keywords{num_keys}='ENFORCE_VERG';			    ENFORCE_VERG=num_keys; 			num_keys=num_keys+1;
one_time_keywords{num_keys}='NULL_VALUE';				NULL_VALUE=num_keys; 			num_keys=num_keys+1;
one_time_keywords{num_keys}='PATCH_OFF';				PATCH_OFF=num_keys;	 			num_keys=num_keys+1;
one_time_keywords{num_keys}='AD_RANGE';				    AD_RANGE=num_keys; 				num_keys=num_keys+1;
one_time_keywords{num_keys}='X_DEG_FULL_SCALE';		    X_DEG_FULL_SCALE=num_keys; 	    num_keys=num_keys+1;
one_time_keywords{num_keys}='Y_DEG_FULL_SCALE';		    Y_DEG_FULL_SCALE=num_keys; 	    num_keys=num_keys+1;
one_time_keywords{num_keys}='FIXED_SEED';				FIXED_SEED=num_keys; 			num_keys=num_keys+1;
% added additional params for screen and interocular distance MLM 7/2/2002
one_time_keywords{num_keys}='SCREEN_XSIZ';              SCREEN_XSIZ=num_keys;           num_keys=num_keys+1;
one_time_keywords{num_keys}='SCREEN_YSIZ';              SCREEN_YSIZ=num_keys;           num_keys=num_keys+1;
one_time_keywords{num_keys}='VIEW_DIST';                VIEW_DIST=num_keys;             num_keys=num_keys+1;
one_time_keywords{num_keys}='IO_DIST';                  IO_DIST=num_keys;               num_keys=num_keys+1;
one_time_keywords{num_keys}='DOT_SIZE_GL';              DOT_SIZE_GL=num_keys;           num_keys=num_keys+1;
one_time_keywords{num_keys}='SOFTWARE_CALIB_STATUS';    SOFTWARE_CALIB_STATUS=num_keys; num_keys=num_keys+1;
%%%% added subject type BJP 1/3/01
one_time_keywords{num_keys}='SUBJECT_TYPE';			    SUBJECT_TYPE=num_keys;		    num_keys=num_keys+1;
%Note: Keep PROTOCOL down here at bottom so all keys are allocated into array.  This is for backwards compatibility.
one_time_keywords{num_keys}='PROTOCOL';				    PROTOCOL=num_keys; 				num_keys=num_keys+1;

NUM_ONE_TIME_PARAMS = num_keys - 1;


%A list of keywords for eye calibration parameters that only get dropped after the FIRST trial
num_keys = 1;
eye_calib_keywords{num_keys}='SOFT_CAL_LEYE_HORIZ';	 SOFT_CAL_LEYE_HORIZ=num_keys; 				num_keys=num_keys+1;
eye_calib_keywords{num_keys}='SOFT_CAL_LEYE_VERT';	    SOFT_CAL_LEYE_VERT=num_keys; 				num_keys=num_keys+1;
eye_calib_keywords{num_keys}='SOFT_CAL_REYE_HORIZ';	 SOFT_CAL_REYE_HORIZ=num_keys; 				num_keys=num_keys+1;
eye_calib_keywords{num_keys}='SOFT_CAL_REYE_VERT';	    SOFT_CAL_REYE_VERT=num_keys; 				num_keys=num_keys+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% added binding stuff BJP 1/3/01

%obj parameters
num_keys = 1;
obj_keywords{num_keys} = 'OBJ_STATUS';					OBJ_STATUS=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_DIR';						OBJ_DIR=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_PHASE';					OBJ_PHASE=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_SKEW_AMPL';				OBJ_SKEW_AMPL=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_SKEW_ANGLE';			    OBJ_SKEW_ANGLE=num_keys;		num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_TRANS_X';				    OBJ_TRANS_X=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_TRANS_Y';				    OBJ_TRANS_Y=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_ROTATION';				OBJ_ROTATION=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_SCALE_X';				    OBJ_SCALE_X=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_SCALE_Y';				    OBJ_SCALE_Y=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_HDISP';					OBJ_HDISP=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_JITTER';					OBJ_JITTER=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_LUM_MULT';				OBJ_LUM_MULT=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_EDGE_WIDTH';			    OBJ_EDGE_WIDTH=num_keys;		num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_AMPL';					OBJ_AMPL=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_AP_TYPE';				    OBJ_AP_TYPE=num_keys;			num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_START_DELAY';			    OBJ_START_DELAY=num_keys;		num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_FREQ';					OBJ_FREQ=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_AMPL_RATIO';			    OBJ_AMPL_RATIO=num_keys;		num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_RLUM';					OBJ_RLUM=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_GLUM';					OBJ_GLUM=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_BLUM';					OBJ_BLUM=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_START_OFFSET';			OBJ_START_OFFSET=num_keys;		num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_STOP_OFFSET';			    OBJ_STOP_OFFSET=num_keys;		num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_TRAJ_ORIENT';			    OBJ_TRAJ_ORIENT=num_keys;		num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_SPEED';					OBJ_SPEED=num_keys;				num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_TRANS_AXIS';				OBJ_TRANS_AXIS=num_keys;		num_keys=num_keys+1;
obj_keywords{num_keys} = 'OBJ_OPEN_PROB';				OBJ_OPEN_PROB=num_keys;		    num_keys=num_keys+1;

NUM_OBJ_PARAMS = num_keys - 1; %use last number

%BAR PARAMETERS
num_keys = 1;
bar_keywords{num_keys} = 'BAR_STATUS';					BAR_STATUS=num_keys;			num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_DIR';						BAR_DIR=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_FREQ';					BAR_FREQ=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_AMPL';					BAR_AMPL=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_X_CTR';					BAR_X_CTR=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_Y_CTR';					BAR_Y_CTR=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_LENGTH';					BAR_LENGTH=num_keys;			num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_TYPE';					BAR_TYPE=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_WIDTH';					BAR_WIDTH=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_ANGLE';					BAR_ANGLE=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_RLUM';					BAR_RLUM=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_GLUM';					BAR_GLUM=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_BLUM';					BAR_BLUM=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_PHASE';					BAR_PHASE=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_OCCLA';					BAR_OCCLA=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_OCCLB';					BAR_OCCLB=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_HDISP';					BAR_HDISP=num_keys;				num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_OCCLWID';				    BAR_OCCLWID=num_keys;			num_keys=num_keys+1;
bar_keywords{num_keys} = 'BAR_OCCLANG';				    BAR_OCCLANG=num_keys;			num_keys=num_keys+1;

NUM_BAR_PARAMS = num_keys - 1; %use last number

%bkgnd parameters
num_keys = 1;
bkgnd_keywords{num_keys} = 'BKGND_BACK_COLOR';		BKGND_BACK_COLOR=num_keys;		num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_DOT_COLOR';		BKGND_DOT_COLOR=num_keys;		num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_AP_RLUM';			BKGND_AP_RLUM=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_AP_GLUM';			BKGND_AP_GLUM=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_AP_BLUM';			BKGND_AP_BLUM=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_BIN_CORR';		BKGND_BIN_CORR=num_keys;		num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_FLICKER';			BKGND_FLICKER=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_DWELL';			BKGND_DWELL=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_HDISP';			BKGND_HDISP=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_DISP_TYPE';		BKGND_DISP_TYPE=num_keys;		num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_DOT_DENSITY';		BKGND_DOT_DENSITY=num_keys;	    num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_XCTR';			BKGND_XCTR=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_YCTR';			BKGND_YCTR=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_XSIZE';			BKGND_XSIZE=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_YSIZE';			BKGND_YSIZE=num_keys;			num_keys=num_keys+1;
bkgnd_keywords{num_keys} = 'BKGND_DOT_SIZE';		BKGND_DOT_SIZE=num_keys;		num_keys=num_keys+1;

NUM_BKGND_PARAMS = num_keys - 1; %use last number

%below, we define a list of analysis options for each experimental Protocol.
%You can easily add the options for new analyses here...  GCD, 1/2/2000

%%%%DOTS PROTOCOLS%%%%
analysis_strings{DIRECTION_TUNING + 1} = ...
    {	'Plot Tuning Curve',
    'Plot Event Times',
    'Plot Spike Rasters',
    'Plot PSTH',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot ISI Histogram',
    'Plot Joint PSTH',
    'Plot Response Distributions',   
    'Experiment Playback'
};
analysis_strings{SPEED_TUNING + 1} = ...
    {	'Plot Tuning Curve',
    'Plot Histograms',
    'Plot Event Times',
    'Plot PSTH',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot ISI Histogram',
    'Plot Joint PSTH',
    'Experiment Playback',
    'Fit LogGauss (Harris)',
    'Fit Gaussian'
};
analysis_strings{HDISP_TUNING + 1} = ...
    {	'Plot Tuning Curve',
    'Fit Gabor/Sine/Gaussian',
    'Gabor Phase/Position Analysis',
    'Experiment Playback',
    'Gabor Fit Multi-Speeds',
    'Gabor Fit Multi-Speed Separable',
    'Plot Time Evolution of Tuning',
    'Multi-Speed Autocorrelation Analysis',
    'Plot Spike Rasters',
    'Plot Vergence Data',
    'Plot Event Times',
    'Plot PSTH',
    'Compile Stationary/Moving PSTH'
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot ISI Histogram',
    'Plot Joint PSTH',
    'Plot Eye Traces',
    'Plot STSA'
    'Signal/Noise Correlation'
};
analysis_strings{DEPTH_DISCRIM + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Psychometric Only',
    'Response Distributions',
    'Compute Choice Probabilities',
    'Vergence Choice Prob',
    'Time Block Analysis',
    'Experiment Playback',
    'Plot Rasters/Histograms',
    'Plot Microstim',
    'Variance Analysis',
    'Poisson Spike Simulation',
    'Noise Correlation'  ,
    'Plot Cross Correlograms',
    'Behavioral Bootstrap'   };
analysis_strings{DIRECTION_DISCRIM + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histogram',
    'Behavioral Bootstrap'  };
analysis_strings{HDISP_GRAD_TUNING + 1} = ...
    {	'Plot Grad Angle Tuning Curve',
    'Simulate Gradient Response',
    'Fit with Sine Wave',
    'Seq F-Test for Fixed/Unfixed Freq in Sin Fit', 
    'Fit with Multiple Sine Waves',
    'Simultaneous Fit with Phase Slop', 
    'Simulate Heterogeneous Surround',
    'Sorted Spike Correlogram',
    'Plot Slant Tuning Curve',
    'Plot Vergence Data'	};
analysis_strings{SIZE_TUNING + 1} = ...
    {	'Plot Tuning Curve',
    'Plot Rasters/Histograms',
    'Plot Vergence Data',
    'Plot Auto and Cross Correlograms'};
analysis_strings{SURROUND_MAPPING + 1} = ...
    {	'Plot Surround Response Map',
    'Plot Rasters/Histograms',
    'Plot Vergence Data'};
analysis_strings{SURROUND_TUNING + 1} = ...
    {	'Plot Surround Tuning Curve',
    'Plot Rasters/Histograms',
    'Plot Vergence Data'};
analysis_strings{RELATIVE_DISPARITY + 1} = ...
    {	'Plot Tuning Curves',
    'Plot Vergence Angle',
    'Plot Rasters/Histograms',
    'Fit Gabor'};
analysis_strings{STEREOACUITY + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Response Distributions',
    'Compute Choice Probabilities',
    'Compute Choice Probabilities from Vergence',
    'Plot Rasters/Histograms',
    'Plot Microstim',
    'Behavioral Bootstrap'  };
analysis_strings{EYE_CALIBRATION + 1} = ...
    {	'Linear Calibration',
    'Nonlinear Calibration'
    'Compare Calibrations'};
analysis_strings{RF_MAPPING + 1} = ...
    {	'Fit with 2D Gaussian' };
analysis_strings{TRANS_RELATIVE_DISPARITY + 1} = ...
    {	'Plot Tuning Curves',
    'Plot Vergence Angle',
    'Plot Rasters/Histograms',
    'Fit Gabor'};
analysis_strings{SIM_DIST_DISP_ONLY + 1} = ...
    {	'Plot Tuning Curves',
    'Plot Rasters/Histograms',
    'Fit Gabor'};
analysis_strings{DEPTH_DISCRIM_NOVAR + 1} = ...
    {	'Plot Psychometric Var/Novar',
    'Plot Rasters/Histograms',
    'Poisson Spike Simulation' };
analysis_strings{SIM_DIST_VERG_ONLY + 1} = ...
    {	'Plot Vergence Data'};
analysis_strings{SIM_DIST_DISP_VERG + 1} = ...
    {	'Plot Tuning Curves',
    'Plot Rasters/Histograms',
    'Fit Gabor',
    'Plot Vergence Data',
    'Refit Eye Data'};
analysis_strings{SIM_DIST_CURVATURE_DISCRIM + 1} = ...
    {	'Plot Psychometric'};
analysis_strings{AXIS_CUED_DIRECTION_TUNING + 1} = ...
    {   'Plot Psychometric'};
analysis_strings{CUED_DIRECTION_DISCRIM + 1} = ...
    {   'Plot Psychometric',
    'Plot 1d Psychometric',
    'Bootstrap Curve Fits',
    'Plot Neurometric/Psychometric',
    'Cue Direction Effects',
    'PSTH',
    'Compute CP'};
analysis_strings{ORIENT_CUE_DIRECTION_TUNING + 1} = ...
    {   'Plot Behavior'};

%---------------------------------------------------------------
%Protocol analysis strings for MOOG Protocols, added 6/27/03 GCD
%---------------------------------------------------------------
analysis_strings{DIRECTION_TUNING_3D + 1} = ...
    {	'Plot Tuning Surface',
    'Fit Optic Flow Tuning (Zack)',
    'Plot Event Times',
    'Plot Spike Rasters',
    'Plot PSTH',
    'Plot PSTH_Anuk',
    'Run PSTH_Gaussfit',
    'Eye trace',    
    'Plot Tuning Azimuth',
    'MU activity'
};
analysis_strings{HEADING_DISCRIM + 1} = ...
    {   'Plot Psycho_neuro_cum',
        'Plot Psycho_neurometric',        
        'Plot Psychometric',
        'Plot Accelerometer',
        'Plot Accelerometer_cum',
        'Plot Response Distribution',
        'Plot CP Distribution',
        'Plot CP Accelerometer',
        'Plot bootstp',
};
analysis_strings{HEADING_DISCRIM_FIXONLY + 1} = ...
    {   'Plot Heading_tuning',
};
analysis_strings{DIR3D_VESTIB_ACCEL + 1} = ...
    {	'Plot Tuning Surface_ves',
    'Plot PSTH',
};
analysis_strings{DIR3D_VARY_FIXATION + 1} = ...
    {	'Plot Fixation Tuning Surface',
    'Plot Fixation PSTH',
    'Output Data for Curve Fitting',
    'Fit Optic Flow Tuning (Zack)',
    'Plot PSTH_Anuk_fix'
};
analysis_strings{AZ_1D_VARY_FIXATION + 1} = ...
    {	'Plot 1D Fixation Tuning Curves_Azimuth',
    'Plot 1D Fixation PSTH'
};
analysis_strings{EL_1D_VARY_FIXATION + 1} = ...
    {	'Plot 1D Fixation Tuning Curves_Elevation',
    'Plot 1D Fixation PSTH'
};
analysis_strings{DIR2D_CUE_CONFLICT + 1} = ...
{	
	'Plot Tuning Surface_conflict',
    'Fit_cue_conflict_2D'
};
analysis_strings{DIRECTION_REVCORR + 1} = ...
{
    'Direction Reverse Correlation',
    'Mapping Receptive Field'
    'Direction Tuning'
};
analysis_strings{AZIMUTH_TUNING_1D + 1} = ...
{
    'Plot Tuning Azimuth'
};
analysis_strings{ROTATION_TUNING_3D + 1} = ...
{
    'Plot Rotation Tuning 3D',
    'Plot Rotation PSTH',
    'Plot Eye Trace',
    'Plot PSTH_Anuk_rotation',
};
analysis_strings{RVOR_PURSUIT + 1} = ...
{
    'Plot Eye Trace',
    'Plot Pursuit Direction Tuning'
};

analysis_strings{TILT_TRANSLATION + 1} = ...
{
    'Plot  TILT_TRANSLATION PSTH',
    'Plot Direction tuning'
};%added by AHC 
%-------------------------------------------
%Protocol Strings for Surface Orientation Stuff - JDN  01/22/04
%-------------------------------------------
analysis_strings{SURF_DEPTH_TUNING+1} = ...
    {'Plot Tuning Curve'};
analysis_strings{SURF_SPEED_TUNING+1} = ...
    {'Plot Tuning Curve'};
analysis_strings{SURF_RF_MAPPING+1} = ...
    {'Plot Tuning Curve'};
analysis_strings{SURF_DIRECTION_TUNING+1} = ...
    {'Plot Tuning Curve'};
analysis_strings{SURF_SIZE_TUNING+1} = ...
    {'Plot Tuning Curve'};
analysis_strings{SURF_TUNING+1} = ...
    {'Plot Tuning Curve'
     'Fit with Ind. Sin Fxn'
     'Plot Vector Average'
     'Plot Vector Average normalized by STD'
     'Plot Vector Average with Tuning Bias (STD)'};
 analysis_strings{SLANT_TUNING+1} = ...
     {'Plot Tilt Slant Tuning Curves'};
     

%-------------------------------------------
%Protocol Strings for Motion Parallax - JWN  09/16/04
%-------------------------------------------
analysis_strings{MOTION_PARALLAX_FIX+1} = ...
{   'Plot Depth-from-MP Tuning Curves',
    'Plot Depth-from-MP Tuning Curves EO/HO', 
    'Plot PSTHs',
    'Get Modulation Indices',
    'Plot Eye Movement Data',
    'Plot Eye Response Correlations'
    'Analyze RF Location'
};


%%%%FOR GRATINGS%%%%

analysis_strings{GRATING_ORIENTATION+1} = ...
{	'Plot Tuning Curve',
    'Plot Spike Rasters',
    'Plot PSTH',
    'Pooling'
};

analysis_strings{GRATING_SPATIAL_FREQ+1} = ...
{	'Plot Tuning Curve',
    'Plot Spike Rasters',
    'Plot PSTH',
    'PSTH Fourier Analysis',
    'Stimulus Onset',
    'Pooling'
};

analysis_strings{GRATING_TEMPORAL_FREQ+1} = ...
{	'Plot Tuning Curve',
    'Plot Spike Rasters',
    'Plot PSTH',
    'Pooling'
};

analysis_strings{GRATING_CONTRAST+1} = ...
{	'Plot Tuning Curve',
    'Plot Spike Rasters',
    'Plot PSTH',
    'Pooling'
};

analysis_strings{GRATING_SIZE+1} = ...
{	'Plot Tuning Curve',
    'Plot Spike Rasters',
    'Plot PSTH',
    'Pooling'
};

analysis_strings{MOTION_COHERENCE+1} = ...
{	'Plot Tuning Curve',
    'Plot Spike Rasters',
    'Plot PSTH'
};

analysis_strings{GRATING_RF_MAP+1} = ...
{	'Fit 2D Gaussian',
    'Plot Spike Rasters',
    'Plot PSTH'
};


%% ADDED BINDING PROTOCOL STRINGS BJP 1/3/01
analysis_strings{FIXATION + 1} = ...
{	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',  
    'Plot Auto and Cross Correlograms'
};

analysis_strings{FIX_1_23_45 + 1} = ...
{	'Plot Neurometric/Psychometric',
    'Plot Event Times',
    'Plot PSTH',
    'Compute Choice Probability',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Compare Correlation Levels',
    'Plot ISI Histogram',
    'Plot Joint PSTH',
    'Plot Response Distributions',
    'Plot Eye Traces',
    'Plot Rasters',
    'Plot LFPs',
    'Plot Spike-Triggered Averages',
    'Plot Correlogram Timecourse',
    'Plot Coherence Spectrogram',
    'Simulate Spike Trains'
};

analysis_strings{FIX_1_23 + 1} = ...
{	'Plot Neurometric/Psychometric',
    'Plot Event Times',
    'Plot Rasters'
    'Plot PSTH',
    'Compute Choice Probability',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Compare Correlation Levels',
    'Plot ISI Histogram',
    'Plot Joint PSTH',
    'Plot Response Distributions',
    'Plot Eye Traces',
    'Plot LFPs',
    'Plot Spike-Triggered Averages',
    'Plot Correlogram Timecourse',
    'Plot Coherence Spectrogram',
    'Simulate Spike Trains'
};

analysis_strings{BIND_DIR_TUNING + 1} = ...
    {	'Plot Direction Tuning Curve',
    'Plot Rasters',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs',
    'Plot PSTH'
    
};

analysis_strings{BIND_TEMPFREQ_TUNING + 1} = ...
    {	'Plot Temporal Frequency Tuning Curve',
    'Plot Rasters',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs',
    'Plot PSTH'
};

analysis_strings{BIND_SPATFREQ_TUNING + 1} = ...
    {	'Plot Spatial Frequency Tuning Curve',
    'Plot Rasters',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs',    
    'Plot PSTH'
    
};

analysis_strings{BIND_RF_MAPPING + 1} = ...
    {	'Fit with 2D Gaussian',
    'Plot Rasters',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs',    
    'Plot PSTH'
};

analysis_strings{BIND_HDISP_TUNING + 1} = ...
    {	'Plot H Disp Tuning Curve',
    'Plot Rasters',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs',    
    'Plot PSTH'
};

analysis_strings{BIND_BAR_DIR_TUNING + 1} = ...
    {	'Plot Direction Tuning Curve',
    'Plot Rasters',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs',
    'Plot PSTH',
    'Plot Spike-Triggered Averages'
};

analysis_strings{BIND_BAR_SPEED_TUNING + 1} = ...
    {	'Plot Speed Tuning Curve',
    'Plot Rasters',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs',
    'Plot PSTH'
    
};
analysis_strings{FIX_VARY_BACKCOLOR + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Event Times',
    'Plot PSTH',
    'Compute Choice Probability',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Compare Correlation Levels',
    'Plot ISI Histogram',
    'Plot Joint PSTH',
    'Plot Response Distributions',
    'Plot Eye Traces',
    'Plot Rasters',
    'Plot LFPs',
    'Plot Spike-Triggered Averages',
    'Plot Correlogram Timecourse',
    'Plot Coherence Spectrogram',
    'Simulate Spike Trains'
};

analysis_strings{FIX_VARY_HISTORY + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Event Times',
    'Plot PSTH',
    'Compute Choice Probability',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Compare Correlation Levels',
    'Plot ISI Histogram',
    'Plot Joint PSTH',
    'Plot Response Distributions',
    'Plot Eye Traces',
    'Plot Rasters',
    'Plot LFPs',
    'Plot Spike-Triggered Averages',
    'Plot Correlogram Timecourse',
    'Plot Coherence Spectrogram',
    'Simulate Spike Trains'
};


analysis_strings{BIND_DIR_DISC + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',  
    'Plot Auto and Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs'
};

analysis_strings{BIND_BACKCOLOR + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',
    'Plot Auto and Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs'
};
analysis_strings{BIND_BACKHDISP + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',
    'Plot Auto and Cross Correlograms'};
analysis_strings{BIND_OBJ_LUM + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',
    'Plot Auto and Cross Correlograms'};
analysis_strings{BIND_OBJ_WIDTH + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',
    'Plot Auto and Cross Correlograms'};

analysis_strings{BIND_AP_RLUM + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',
    'Plot Auto and Cross Correlograms'};

analysis_strings{BIND_OBJ_POSITION + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',
    'Plot Auto and Cross Correlograms'};
analysis_strings{BIND_DOT_DENSITY + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',
    'Plot Auto and Cross Correlograms'};

analysis_strings{BIND_COHER_DISC + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',  
    'Plot Auto and Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs'
};

analysis_strings{BIND_CLOSURE_DISC + 1} = ...
    {	'Plot Neurometric/Psychometric',
    'Plot Microstim Effect',
    'Compute Choice Probability',
    'Plot Rasters/Histograms',  
    'Plot Auto and Cross Correlograms',
    'Plot Eye Traces',
    'Plot LFPs'
};

analysis_strings{FIX_TRANS_QUAD + 1} = ...
    {	'Plot Speed Tuning Curve',
    'Plot Rasters',
    'Plot Auto Correlograms',
    'Plot Cross Correlograms',
    'Plot Eye Traces',
    'Plot PSTH',
    'Plot Circ Trans Tuning Curve'
};