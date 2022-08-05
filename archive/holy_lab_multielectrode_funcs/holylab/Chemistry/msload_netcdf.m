function ncdata = msload_netcdf(filename)
% MSLOAD_NETCDF: load mass spec .cdf files from Gross lab
%
% The core data are stored as a long vector, indexed as "points". Each
% "scan" consists of a sequence of many points.
%
% Syntax:
%   ncdata = msload_netcdf(filename)
% where
%   filename is a string containing the name of the .cdf file
% and
%   ncdata is a structure with mostly self-explanatory fields. Here are the
%   ones that might need a bit of explanation:
%     scan_time: the start time of each scan, in seconds, relative to file
%       start
%     scan_index: scan_index(i) is the point number at the beginning of the
%       ith scan
%     scan_mzmin/scan_mzmax: the min and max mz for each scan
%     mz: arranged as a vector, the mz value at each point
%     intensity: arranged as a vector, the ion count at each point

% Copyright 2009 by Timothy E. Holy

  [ncid,status] = mexnc('OPEN',filename);
  if (status ~= 0)
    error(mexnc('STRERROR',status));
  end
  
  % Get header info
  ncdata.date_time = get_str_attribute(ncid,'experiment_date_time_stamp');
  ncdata.type = get_str_attribute(ncid,'experiment_type');
  ncdata.sample = get_str_attribute(ncid,'sample_comments');
  %sample_prep_comments?
  varid = get_varid(ncid,'scan_acquisition_time');
  [ncdata.scan_time,status] = mexnc('GET_VAR_DOUBLE',ncid,varid);
  varid = get_varid(ncid,'scan_type');
  [ncdata.ms_level,status] = mexnc('GET_VAR_INT',ncid,varid);
  varid = get_varid(ncid,'scan_index');
  [ncdata.scan_index,status] = mexnc('GET_VAR_INT',ncid,varid);
  varid = get_varid(ncid, 'mass_range_min');
  [ncdata.scan_mzmin,status] = mexnc('GET_VAR_DOUBLE',ncid,varid);
  varid = get_varid(ncid, 'mass_range_max');
  [ncdata.scan_mzmax,status] = mexnc('GET_VAR_DOUBLE',ncid,varid);
%   varid = get_varid(ncid, 'resolution');
%   [ncdata.resolution,status] = mexnc('GET_VAR_DOUBLE',ncid,varid);
  
  % Get raw data
  varid = get_varid(ncid,'total_intensity');
  [ncdata.total_intensity,status] = mexnc('GET_VAR_DOUBLE',ncid,varid);
  varid = get_varid(ncid,'mass_values');
  [ncdata.mz,status] = mexnc('GET_VAR_FLOAT',ncid,varid);
  varid = get_varid(ncid,'intensity_values');
  [ncdata.intensity,status] = mexnc('GET_VAR_INT',ncid,varid);
  
  mexnc('CLOSE', ncid);
end

function varid = get_varid(ncid,varname)
  [varid,status] = mexnc('INQ_VARID',ncid,varname);
  if (status ~= 0)
    error(mexnc('STRERROR',status));
  end
end

function str = get_str_attribute(ncid,attrname)
  [str,status] = mexnc('GET_ATT_TEXT',ncid,-1,attrname);
  if (status ~= 0)
    error(mexnc('STRERROR',status));
  end
  if (str(end) == 0)
    str = str(1:end-1);
  end
  str = strtrim(str);
end
