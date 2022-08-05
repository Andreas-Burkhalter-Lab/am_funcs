% General utility functions
%
% Finding and indexing:
%   findainb          - Return index vector for elements of a in b.
%   indexainb         - Detemine which members of a are in b.
%   time2indx         - Version of findainb with approximation.
%   agglabel          - Aggregate integer labels efficiently.
%   anyisnan          - Check whether any entry in coordinate slices is NaN.
%   countnans         - Count the number of NaNs in coordinate slices
%   count2index       - Calculate the indices given a series counts
%   find_first_slot   - Find the first empty slot in a sequence
%   increment_counter - Advance a multi-digit counter (with carry)
%   ind2sub_matrix    - Like ind2sub, but coordinates are specified as a matrix
%   sub2ind_matrix    - Like sub2ind, but coordinates are specified as a matrix
%   make_counting_index - Create a vector to index variable-length blocks
%
% String manipulation:
%   exclude_strings  - eliminate strings from a cell array that match pattern
%   key2value        - retrieve the value by key in a string that includes many key=value pairs
%   sectionkey2value - retrieve the value by section&key from a string in INI file format
%   update_value     - update the value by key in a string that includes many key=value pairs
%   split_str        - split a string (a label list) into string list
%   split_label      - split the label string(from a Merec file header) to labels
%   segment          - similar to split_str(), this func returns index ranges 
%   segment_indices  - simplify discrete indices into index ranges
%   split_dbl        - split a string into double list
%   split_int        - split a string into int list
%   string_endswith  - test if a string ends w/ a pattern
%   concstring2molar - convert string concentrations to numerical units
%   safe_eval        - eval a string w/o worrying about corrupting current namespace
%
% Array manipulation:
%   squeeze2       - Remove singleton dimensions, even from 2-D arrays
%   coords2spans      - convert coordinate lists into ranges
%   spans2coords      - convert ranges into coordinate lists
%   array_snip_common - compute coordinates for snipping out regions of arrays
%
% Structure manipulation:
%   copyfields     - Copy field values from one structure to another
%   extract_fields - Convert structures into vectors
%   fill_fields    - Convert vectors into structures
%
% Function handles:
%   foreach_g      - Apply a function handle to an array of arguments
%
% Graphics handles:
%   free           - Safely delete handles (won't complain for non-handles)
%   isvisible      - test if the GUI objects are visible
%
% User interaction:
%   keystepper        - use keypresses to cycle through values
%   progress          - create/update/destroy a progress bar
%   progress_bar      - create/update/destroy a progress bar
%   nonblock_inputdlg - non-block Input dialog box.
%
% Date and time functions:
%   datenum_g    - a wrapper to datenum to parse more formats
%   extract_time - parse a string for time
%   time_in_secs - parse time strings, return in units of seconds
%
% Memory:
%   sizeof      - Return the size of a builtin data type, in bytes
%   cleard      - Clear all variables except string vars beginning with 'd'
%
% Path manipulation:
%   clean_path    - exclude inappropriate directories from the path
%   make_abs_path - generate absolute path from a path
% 
% File name manipulation:
%   dirbyname      - Select files matching pattern, sorted alphabetically.
%   dirtime        - Select files matching pattern, sorted by creation time.
%   UIGetFiles     - GUI for selecting multiple files.
%   canonicalize_file_name - create full-path version of a file name.
%   find_first_pattern - find the first file pattern that has existent file(s).
%   replace_extension  - replace the filename's extension name
%   replace_filename   - replace the filename's name part (retain the same parent dir)
%   replace_parent_dir - replace the filename's directory part
%
% File information and manipulation:
%   fileexist - test if a file exists
%   filesize  - get the file's size
%   fopen_g   - a wrapper function to fopen() and test if the fid is reused.
%   ismat         - test whether a file is a MATLAB file (.mat)
%   is_merec_file - test if a file was saved by Merec
%   is_new_merec  - test if a file was saved by a "modern" Merec
%   reorder_merec - reorder channels in a .merec file
%   sort_ai_by_time - order .merec, .ssnp, etc files by recording time
%   sortfilebytime  - sort files by modification time (from the earliest to the latest)
%
% File reading and writing:
%   dlmwrite0    - Write "spreadsheet file" like dlmwrite, except handling 0s.
%   load_text_file - return a text file's content as a string
%   MATToText    - Convert simple .mat files to ASCII
%   readheader   - Read the header of all types of electrophysiology data files
%   read_ascii_header    - read the Ascii header of Merec data, envelope data, etc
%   read_envelope_header - read envelope header into a matlab structure
%   read_merec_header    - read Merec header into a matlab structure
%   supportlegacyheaders - a flag determining data format compatibility
%   resolvelegacyheader  - a gateway which deals with old header types
%   ReadAIHeaderHarvard  - Read the header from legacy raw waveform files
%   WriteAIHeaderHarvard - Write the header for legacy raw waveform files
%   READAIHEADERWU1      - Read the header of raw waveform files (pre-merec)
%   WRITEAIHEADERWU1     - Write the header for legacy raw waveform files (pre-merec)
%   update_header        - update header content
%   update_magic         - update magic sequence that tells file type
%   slice_merec          - copy parts of some .merec files and form a new .merec file
%
% MEX files:
%   mexpp                - Simplifies building sophisticated C++ MEX files
%   MatlabTraits.h       - Traits classes that simplify MEX creation
%   MatlabIO.h           - Variable conversion (matlab/C/STL/Eigen) and .mat file read/write
%   BoostThreadPool.h    - Task-based processing for multithreading, with progress and CTRL-C interrupts
%   BoostThreadPoolExample.cpp   - An example file using BoostThreadPool.h
%   MatlabStubs.h        - Convert common mex* functions into C
%   mat_variables.h      - Older approach to reading/writing from .mat files
%   matlab_arg.h         - Another approach to matlab/C variable conversion
%   threadpool.hpp       - Another approach for task-based processing
%   matlabarg.h          - A convenience wrapper (alternative to MatlabIO/MatlabTraits)
%   timer_g.h            - Measuring timing in MEX files (or use Boost.Time)
%
% Miscellaneous (should probably go somewhere else!):
%   CellToMatrix        - Convert cell array of vectors to a matrix.
%   IntersectIntervals  - Intersection of 2 intervals on real line.
%   sqrdist             - Euclidean distances between pairs of points.
%   sqrdistslow         - Memory-safe alternative of sqrdist.
%   make_vector         - convert the input into a row or column vector
%   unit_tform          - return a unit TForm
%   sleep_g             - make the caller sleep for specified time
%   noneven_mean        - calculate mean values even there're missing values
%   show_help           - show a dialog w/ help message
%
% Large file systems (deprecated):
%   should_use_lfs  - test if we have to use our LFS wrapper
%   openlfs         - open a file with LFS mode
%   closelfs        - close a file opened by openlfs()
%   readcharlfs     - read a block of chars into buf
%   readint16lfs    - read a block of int16's into wave
%   readuint16lfs   - read a block of uint16's into wave
%
% high density array (HDA)
%   get_hda_chan    - return HDA channel numbers (in Comedi's numbering)
%   get_hda_holylab - return HDA Comedi numbers as a matrix organized by their physical locations
%   get_high_density_array - return HDA channel numbers (in LabView's numbering)
%   merec2berrydata - convert .merec recorded in wrong layout to the right layout
%   merec2berrydata_batch - a wrapper to merec2berrydata to make batch processing easier


==> edgereflect.m <==
function xo = edgereflect(sz,x)
% EDGEREFLECT: coordinates of arrays with mirror-reflection symmetry

==> fillmask.m <==
function xf = fillmask(x,isgood)
% FILLMASK: replace masked values by interpolated ones


==> interpmatrix.m <==
function M = interpmatrix(varargin)
% INTERPMATRIX: multidimensional linear interpolation

===> ndgrid_as_matrix.m <==
function gridout = ndgrid_as_matrix(varargin)
% NDGRID_AS_MATRIX: give grid coordinates in column format


==> pwlin_bounds.m <==
function yi = pwlin_bounds(x,y,xi)
% PWLIN_BOUNDS: implement a piecewise linear monotonic function, respecting bounds


==> split_into_contiguous_regions.m <==
function [region,regionIndex] = split_into_contiguous_regions(index,irange)
% SPLIT_INTO_CONTIGUOUS_REGIONS: snippet so that overlapping regions are combined

==> sqrdistslow.m <==
function d2 = sqrdistslow(X)
% SQRDISTSLOW: compute matrix of square distances in a slow but memory efficient way

