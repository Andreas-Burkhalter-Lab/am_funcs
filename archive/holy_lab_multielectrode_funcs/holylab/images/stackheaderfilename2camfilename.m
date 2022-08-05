function camfile=stackheaderfilename2camfilename(headerfile, parsedheader)
if parsedheader.header_version <= 3.0
    camfile = [headerfile '.dat'];
    if(fileexist(camfile))
        return
    end
    camfile=replace_extension(headerfile, '.dat');
else
    camfile=[headerfile '.cam'];
    if(fileexist(camfile))
        return
    end
    camfile=replace_extension(headerfile, '.cam');
end

