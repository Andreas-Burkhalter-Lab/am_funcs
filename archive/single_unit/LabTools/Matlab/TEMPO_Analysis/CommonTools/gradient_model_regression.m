function gradient_model_regression

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*'};
lines = {'b-' 'r-' 'g-' 'k-' 'm-'};
lines2 = {'b--' 'r--' 'g--' 'k--' 'm--'};

PATH = 'Z:\Users\Jerry\GradAnalysis\';
BATCH = 'measured_model_batch.m';
BATCHIN = [PATH BATCH];

batchid = fopen(BATCHIN);

while 1
    FILE = fgetl(batchid);
    if ~ischar(FILE), break, end

    FILEIN = [PATH FILE];
    f = strfind(FILE, 'for');
    datafilename = FILE(f+4:length(FILE));
    
    fid = fopen(FILEIN);
    
    i = 1;
    while 1
        datain{i} = fgetl(fid);
        if ~ischar(datain{i}), break, end
        i = i+1;
    end
    
    measured_and_model_data = zeros(i-1, 5);
    
    for j = 1:i-1
        measured_and_model_data(j,:) = str2num(datain{j});
    end
    
    mean_shifted_data = measured_and_model_data;
    mean_measured = mean(measured_and_model_data(:,3));
    mean_model = mean(measured_and_model_data(:,5));
    
    mdisp = munique(measured_and_model_data(:,1));
    
    outfile = [PATH 'r_stats_for_gradient_model_each_mdisp.dat' ];
    fid2 = fopen(outfile, 'a');
    
    display = 0;
    if display
        temp_measured = figure;
        temp_model = figure;
    end
    rlist = zeros(length(mdisp), 1);
    for i=1:length(mdisp)
        mdisp_ind = find(measured_and_model_data(:,1)==mdisp(i));
        one_list = ones(length(mdisp_ind), 1);
        X = [one_list measured_and_model_data(mdisp_ind, 5)];
        [b,bint,r,rint,stats] = regress(measured_and_model_data(mdisp_ind, 3),X);
        rlist(i) = stats(1);
        r_val(i) = sqrt(rlist(i));
        if b(2) < 0
            r_val(i) = -r_val(i);
        elseif b >= 0
            r_val(i) = r_val(i);
        end
        fprintf(fid2, '%s\t%1.2f\t%1.4f\t%1.4f\t%3.4f\t%1.5f\n', datafilename, mdisp(i), r_val(i), stats(1), stats(2), stats(3));
        if display    
            figure(temp_measured)
            hold on
            plot(measured_and_model_data(mdisp_ind, 3), lines{i});
            figure(temp_model)
            hold on
            plot(measured_and_model_data(mdisp_ind, 5), lines{i});
        end
        
        %shift means for measured responses
        mdisp_mean = mean(measured_and_model_data(mdisp_ind, 3));
        mean_diff = mean_measured - mdisp_mean;
        mean_shifted_data(mdisp_ind, 3) = mean_shifted_data(mdisp_ind, 3) + mean_diff;
        
        if display
            figure(temp_measured)
            hold on
            plot(mean_shifted_data(mdisp_ind, 3), lines2{i});
        end
        
        %shift means for model responses        
        mdsip_mean = mean(measured_and_model_data(mdisp_ind, 5));
        mean_diff = mean_measured - mdisp_mean;
        mean_shifted_data(mdisp_ind, 5) = mean_shifted_data(mdisp_ind, 5) + mean_diff;
        
        if display
            figure(temp_model)
            hold on
            plot(mean_shifted_data(mdisp_ind, 5), lines2{i});
        end
    end
    outfile2 = [PATH 'r_stats_for_gradient_model_combined.dat'];
    fid3 = fopen(outfile2, 'a');
    one_list = ones(length(mean_shifted_data),1);
    X = [one_list mean_shifted_data(:,5)];
    [b,bint,r,rint,stats] = regress(mean_shifted_data(:, 3),X);
    r_val(i) = sqrt(stats(1));
    if b(2) < 0
        r_val(i) = -r_val(i);
    elseif b >= 0
        r_val(i) = r_val(i);
    end
    fprintf(fid3, '%s\t%1.4f\t%1.4f\t%3.4f\t%1.5f\n', datafilename, r_val(i), stats(1), stats(2), stats(3));
    
    %average r squared values    
    outfile3 = [PATH 'avg_r_stats_for_gradient_model.dat'];
    fid4 = fopen(outfile3, 'a');

    avg_r_squared = mean(rlist);
    avg_r = mean(r_val);
    fprintf(fid4, '%s\t%1.4f\t%3.4f\n', datafilename, avg_r, avg_r_squared);    
    fclose(fid4);
    
    fclose(fid3); %rstats for combined mean disparities
    fclose(fid2);  %rstats each mean disparity
    fclose(fid); %input data file
end

fclose(batchid);


