function [imdata, vol_tune, act_num] = pdocpi_fread(camdata, data_specs, thisActuator)
% PDOCPI_FREAD: read image data acquired during DM tuning
    % 
    % INPUT
    % the acquisition script generates 3 files, (2 *.cam files and one data specs (.mat) file
    % use one of the .cam file at a time with the data specs file
    %acquisition of phase diversity images on computer extreme
    %
    % OUTPUT
    % imdata will a 3-D or 4-D matrix, 2-D image, 3rd dimension is applied voltage, 4th
    % dimension (if applicable) is actuator number
    % vol_tune is the applied voltage; a single number if accompanied by
    % actuator data; or a n by 52 matrix (consisting of each actuator's applied
    % voltage)
    %
    % Example:
    %
    % (for DM calibration data)
    % [im1data, vol_tune, act_num] = pdocpi_fread('cam1.cam','data_specs');
    %
    % (for data acquisition)
    % [im2data, vol] = pdocpi_fread('cam2.cam','data_specs');
    %
    % (just data)
    % imdata = pdocpi_fread('cam1','data_specs');
    %
    %
    % You can also read the data corresponding to a single actuator like
    % this: (This will work only for dm tuning data)
    %   imdata = pdocpi_fread('cam1.cam','data_specs',37);
    % In this example, 37 is the actuator number (not the "37th actuator").
    % You can also use a structure input for data_specs like this:
    %   imdata = pdocpi_fread('cam1.cam',data_specs,37);
    % where data_specs = load('specfilename');
    
    %
    % written by Diwakar, 2009-06-23, and extended by Timothy E. Holy

    if ~isstruct(data_specs)
      % Load from file
      data_specs = load(data_specs);
    end
    
    if isfield(data_specs,'vol_tune')% this means this is the first dm tuning script ...pre version control
        
        
        data_specs = default(data_specs,'pixel_format','*uint16');
        n_actuators = length(data_specs.act_num);
        n_voltages = length(data_specs.vol_tune);
        vol_tune = data_specs.vol_tune;

        [fid,msg] = fopen(camdata);
        if (fid == -1)
          % Can't open file
          error(msg)
        end

        if (nargin < 3)
          % Read all the data
          imdata = fread(fid,data_specs.pixel_format);
          fclose(fid);  % we close first so errors in reshape don't leave dangling fid
          imdata = reshape(imdata, [data_specs.vidRes([2 1]) n_voltages n_actuators]);
          act_num = data_specs.act_num;
        else
          [act_num,act_indx] = intersect(data_specs.act_num,thisActuator);
          if (length(act_num) ~= 1)
            fclose(fid)
            error('Must supply one acquired actuator number');
          end
          stacksz = prod(data_specs.vidRes) * n_voltages;
          status = fseek(fid,(act_indx-1)*stacksz*sizeof(data_specs.pixel_format),'bof');
          if (status < 0)
            fclose(fid);
            error(['File error seeking to position for actuator ' num2str(thisActuator)])
          end
          imdata = fread(fid,stacksz,data_specs.pixel_format);
          fclose(fid);
          imdata = reshape(imdata,[data_specs.vidRes([2 1]) n_voltages]);
        end
    end
    
    if isfield(data_specs,'version') % this is version 1 data acquisition style

        data_specs = default(data_specs,'pixel_fwrite_format','*uint16');
        n_voltages = size(data_specs.apply_vol,1);
        vol_tune = data_specs.apply_vol;
        
        vid1 = data_specs.vid1;
        vidRes = vid1.VideoResolution;

        [fid,msg] = fopen(camdata);
        if (fid == -1)
          % Can't open file
          error(msg)
        end

      % Read all the data
      imdata = fread(fid,data_specs.pixel_fwrite_format);
      fclose(fid);  % we close first so errors in reshape don't leave dangling fid
      imdata = reshape(imdata, [vidRes([2 1]) n_voltages]);

    end
        
    
end