function find_pref_tilt_model

PATH = 'Z:\Users\Jerry\GradAnalysis\';
PATHCURVE = 'Z:\Data\Tempo\Jerry\Analysis\Grad\';
BATCH = 'Model_shift_mean_batch.m';
BATCHIN = [PATH BATCH];

batchid = fopen(BATCHIN);

while 1
    FILE = fgetl(batchid);
    if ~ischar(FILE), break, end

    FILEIN = [PATHCURVE FILE];
    f = strfind(FILE, '.model');
    datafilename = [FILE(1:f-1) '.htb']
    
    fid = fopen(FILEIN);
    
    i = 1;
    datain = {};
    while 1
        datain{i} = fgetl(fid);
        if ~ischar(datain{i}), break, end
        i = i+1;
    end
    
    in_header = 1;
    disp_val = {};
    mdisp = 1;
    done = 2;
    i = 1;
    while done < length(datain)
        temp = str2num(datain{done});
        done = done + 1;
        if ((length(temp) == 8) & (in_header==1))
            disp_val{mdisp}(i,1:2) = temp(1:2);            
            i = i+1;
        elseif length(temp) == 2
            disp_val{mdisp}(i,1:2) = temp(1:2);
            in_header = 0;
            i = i+1;
        elseif ((length(temp) == 8) & (in_header==0))
            in_header = 1;
            mdisp = mdisp + 1;
            i = 1;
            disp_val{mdisp}(i,1:2) = temp(1:2);
        end
    end
    fclose(fid)
    
    fileout = [PATH 'pref_tilt_model_01.21.03.dat'];
    fid2 = fopen(fileout, 'a');
    pref_tilt = zeros(mdisp, 1);
    string_out = '';
    for i=1:mdisp
        [y, ind] = max(disp_val{i}(:,2));
        pref_tilt(i) = disp_val{i}(ind, 1);
         fprintf(fid2, '%s\t%3.2f\n', datafilename, pref_tilt(i));
     end
    fclose(fid2);
end
fclose(batchid)