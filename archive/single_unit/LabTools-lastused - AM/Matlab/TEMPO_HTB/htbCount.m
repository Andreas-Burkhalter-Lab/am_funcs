function n = htbCount(fid)
%htbCount - Get number of databases in a TEMPO HTB file
%
% IN
%   fid             fid returned by htbOpen() or htbOpenw()
%
% OUT
%   n               Number of databases in HTB file
%
% SEE ALSO
%
%   htbOpenw(), htbClose(), htbGetDBCount(), htbReadHeader(), htbReadData()

offset = 0;                             % Offset of first HTB structure in file
n = 0;
while (1)
    status = fseek(fid, offset + 114, 'bof');   % 114=offset to htb.alloc field
    if (status == -1)
        offset = -1;
        return;
    end
        
    [nbytes, count] = fread(fid, 1, 'uint32');  % Read htb.alloc
	 %fprintf('HTBCOUNT - n=%d offset=%d nbytes=%d count%d\n', n, offset, nbytes, count);    
    
    if ((count ~= 1) | (nbytes <= 0))
        offset = -1;
        return;
    end
    
    offset = offset + nbytes;
    n = n + 1;
end;
return;
