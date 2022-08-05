function test_varputg ( ncfile )
% TEST_VARPUTG
%
% This routine tests VARGETG, VARPUTG
%
% Test 1:  test VARPUTG/VARGETG with double precision data
% Test 2:  test VARPUTG/VARGETG with float data, should not be accepted
% Test 3:  
% Test 010:  test writing a short datum to a double precision variable
% Test 100:  VARPUTG with a bad ncid
% Test 101:  VARGETG with a bad ncid
% Test 102:  VARPUTG with a bad varid
% Test 103:  VARPUTG without a missing imap parameter

mexnc ( 'setopts', 0 );

create_testfile ( ncfile );
test_001 ( ncfile );
test_002 ( ncfile );
test_003 ( ncfile );
test_010 ( ncfile );
test_100 ( ncfile );
test_101 ( ncfile );
test_102 ( ncfile );

fprintf ( 1, 'VARPUTG succeeded\n' );
fprintf ( 1, 'VARGETG succeeded\n' );

return




function create_testfile ( ncfile )


%
% ok, first create this baby.
[ncid, status] = mexnc ( 'create', ncfile, nc_clobber_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end



%
% Create the fixed dimension.  
len_x = 100;
[xdimid, status] = mexnc ( 'def_dim', ncid, 'x', len_x );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[z_double_varid, status] = mexnc ( 'def_var', ncid, 'z_double', nc_double, 1, [xdimid] );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[z_float_varid, status] = mexnc ( 'def_var', ncid, 'z_float', nc_float, 1, [xdimid] );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end



[z_short_varid, status] = mexnc ( 'def_var', ncid, 'z_short', nc_short, 1, [xdimid] );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


eps = 0.01;
status = mexnc ( 'put_att_double', ncid, z_short_varid, 'scale_factor', nc_double, 1, eps );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


status = mexnc ( 'put_att_double', ncid, z_short_varid, 'add_offset', nc_double, 1, 0.00 );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[status] = mexnc ( 'enddef', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


%
% CLOSE
status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return















function test_001 ( ncfile );


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_double_varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = rand(25,1);
status = mexnc ( 'VARPUTG', ncid, z_double_varid, [1], [25], [2], [], input_data );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end


[output_data, status] = mexnc ( 'VARGETG', ncid, z_double_varid, [1], [25], [2] );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

output_data = output_data(:);

d = max(abs(output_data-input_data))';
if (any(d))
	error ( 'values written by VARGET do not match what was retrieved by VARPUT\n'  );
	error ( msg );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return








function test_002 ( ncfile );


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[varid, status] = mexnc('INQ_VARID', ncid, 'z_float');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = single(rand(25,1));
try
	status = mexnc ( 'VARPUTG', ncid, varid, [1], [25], [2], [], input_data );
	error ( 'Succeeded when it should have failed' );
end


status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return









function test_010 ( ncfile );


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = int16(rand(25,1)*100);
try
	status = mexnc ( 'VARPUTG', ncid, varid, [1], [25], [2], [], input_data );
	error ( 'succeeded when it should have failed' );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return









function test_100 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_double_varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = rand(25,1);
try
    status = mexnc ( 'VARPUTG', -100, z_double_varid, [1], [25], [2], [], input_data );
	msg = sprintf ( '%s:  %s:  VARPUTG succeeded with a bad ncid\n', mfilename, testid );
	error ( msg );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return







function test_101 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_double_varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

try
    [output_data, status] = mexnc ( 'VARGETG', -100, z_double_varid, [1], [25], [2] );
	msg = sprintf ( '%s:  %s:  VARGET succeeded with a bad ncid\n', mfilename, testid );
	error ( msg );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return








function test_102 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_double_varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = rand(25,1);
try
    status = mexnc ( 'VARPUTG', ncid, -500, [1], [25], [2], [], input_data );
	msg = sprintf ( '%s:  %s:  VARPUTG succeeded with a bad varid\n', mfilename, testid );
	error ( msg );
end


status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return









%
% VARPUTG needs to have a start, count, stride, and then an empty set 
function test_103 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_double_varid, status] = mexnc('INQ_VARID', ncid, 'z_double');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = rand(25,1);
try
    status = mexnc ( 'VARPUTG', ncid, -500, [1], [25], [2], input_data );
	msg = sprintf ( '%s:  %s:  VARPUTG succeeded without a missing imap position \n', mfilename, testid );
	error ( msg );
end


status = mexnc ( 'close', ncid );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

return



function test_003 ( ncfile )

[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[z_short_varid, status] = mexnc('INQ_VARID', ncid, 'z_short');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

[scale_factor, status] = mexnc('GET_ATT_DOUBLE', ncid, z_short_varid, 'scale_factor');
if ( status ~= 0 ), error ( mexnc('strerror',status) ), end

input_data = rand(100,1);
input_data = input_data(1:2:end);
[r,c] = size(input_data);
status = mexnc ( 'VARPUTG', ncid, z_short_varid, [0 0], [r], [2], [], input_data', 1 );
if ( status ~= 0 )
	ncerr_msg = mexnc ( 'strerror', status );
	msg = sprintf ( '%s:  VARPUTG failed, (%s)\n', mfilename, ncerr_msg );
	error ( msg );
end


[output_data, status] = mexnc ( 'VARGETG', ncid, z_short_varid, [0], [r], [2], [], 1 );
if ( status ~= 0 )
	msg = sprintf ( '%s:  ''%s''\n', mfilename,  mexnc ( 'strerror', status ) );
	error ( msg );
end

output_data = output_data(:);

d = max(abs(output_data-input_data))';
ind = find ( d > scale_factor/2 );
if (any(ind))
	msg = sprintf ( 'values written by VARPUTG do not match what was retrieved by VARGETG\n'  );
	error ( msg );
end


fprintf ( 1, 'VARPUTG with scaling (please don''t do this, it''s bad) succeeded\n' );
fprintf ( 1, 'VARGETG with scaling (please don''t do this, it''s bad) succeeded\n' );




