fid = fopen('2008_05_18_220um_10Hz_1repeat_a_reg2_10x10_100sigma_subtr.cam','a');
for idx = 1:2600
%     if length(num2str(idx)) == 1
%         fnum = ['00' num2str(idx)];
%     elseif length(num2str(idx)) == 2
%         fnum = ['0' num2str(idx)];
%     elseif length(num2str(idx)) == 3
%         fnum = num2str(idx);
%     end
    fwrite(fid, uint16(smm_reg_cp(:,:,idx)), 'uint16');
    fprintf('%d..',idx); if mod(idx,20)==0; fprintf('\n'); end;
end

fclose(fid);