%----------------------------------------------------------------------------------------------------------
%-- Path_Defs.m: Path Definitions for various common and protocol-specific analysis routines.  You can
%--	simply modify BASE_PATH to change the base directory for all routines, e.g., when moving code
%--	to a different machine.  GCD, 1/3/2000
%----------------------------------------------------------------------------------------------------------

%Base path specification for protocol-specific analysis routines
BASE_PATH = 'C:\LabTools\Matlab\TEMPO_Analysis\';

%add path for common tools
junk_str = [BASE_PATH 'CommonTools\'];
addpath(junk_str);

%add paths for individual protocol analysis tools
junk_str = [BASE_PATH 'ProtocolSpecific\DirectionTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SpeedTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\DirectionDiscrim'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\SizeTuning'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\RFMapping'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\GratingOrientation'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\GratingTemporalFreq'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\GratingContrast'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\GratingSize'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\GratingRFMap'];
addpath(junk_str);
junk_str = [BASE_PATH 'ProtocolSpecific\Coherence'];
addpath(junk_str);

