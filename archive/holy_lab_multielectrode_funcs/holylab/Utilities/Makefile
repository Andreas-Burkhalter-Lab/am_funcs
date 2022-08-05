#a simple makefile to make LFS related funcs

all: closelfs readint16lfs openlfs test_fid_returned_by_matlab\
     split_str split_int split_dbl readuint16lfs \
     readint8lfs sleep_g

closelfs: closeLFS.cpp
	mex -o $@ closeLFS.cpp

readint16lfs: readint16LFS.cpp
	mex -o $@ readint16LFS.cpp

readuint16lfs: readint16LFS.cpp
	mex -DUINT16SAMPLE -o $@ readint16LFS.cpp

openlfs: openLFS.cpp matlab_arg.cpp matlab_arg.h
	mex -o $@ openLFS.cpp matlab_arg.cpp 

#read_merec_header_cpp: ReadMerecHeader.cpp parsemerecheader.cpp parsemerecheader.h
#	mex -o $@ ReadMerecHeader.cpp   parsemerecheader.cpp

test_fid_returned_by_matlab:  test_fid_returned_by_matlab.cpp
	mex -o $@ test_fid_returned_by_matlab.cpp

split_str: split_str.cpp matlab_arg.cpp matlab_arg.h
	mex split_str.cpp matlab_arg.cpp

split_int: split_int.cpp matlab_arg.cpp matlab_arg.h
	mex split_int.cpp matlab_arg.cpp

split_dbl: split_dbl.cpp matlab_arg.cpp matlab_arg.h
	mex split_dbl.cpp matlab_arg.cpp

readint8lfs: readint8lfs.cpp
	mex readint8lfs.cpp

sleep_g: sleep_g.cpp
	mex sleep_g.cpp

