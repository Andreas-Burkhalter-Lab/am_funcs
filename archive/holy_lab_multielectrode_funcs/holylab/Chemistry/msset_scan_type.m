function ncdata = msset_scan_type(ncdata,mzmin,mzmax)
  ncdata.scan_type = ncdata.scan_mzmin == mzmin & ...
      ncdata.scan_mzmax == mzmax;
end
  