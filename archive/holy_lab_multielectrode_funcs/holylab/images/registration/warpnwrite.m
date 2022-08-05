function warpnwrite(fid, smm_in, u, idx, usave_out)
% warpnwrite will warp and write a .cam file to an already-open file
%
% Syntax: warpnwrite(fid,smm_in,u,idx,usave_out)
%
% Copyright 2010 Julian P Meeks
    stack_size = smm_in.size;
    header = smm_in.header;
    % set up padmat
    padmat = ones(1, stack_size(2), stack_size(3), 'single');
    if ischar(usave_out.bad_stack_file)
        if exist(usave_out.bad_stack_file,'file')
            smm_bad = stackmm(usave_out.bad_stack_file);
        end
    end
    thisu = u;
    found_bad = intersect(idx,usave_out.bad_stacks);
    if isempty(found_bad)
        tempImg = double(smm_in(:,:,:,idx));
    else % modify input stacks to incorporate "repaired" bad stacks
        bad_idx = find(usave_out.bad_stacks == found_bad);
        tempImg = double(smm_bad(:,:,:,bad_idx));
    end
    u_basis = size(usave_out.rmg.image_grid(1).imFixed);
    if sum(abs(diff([stack_size(1:3); u_basis]))) ~= 0
        correct = true;
        if isfield(usave_out.method,'nanpix')
            nanpix = usave_out.method.nanpix;
            nanpad{1} = zeros([stack_size(1) nanpix(1) stack_size(3)]);
            nanpad{2} = zeros([nanpix(2) stack_size(2)+nanpix(1)*2 stack_size(3)]);
            nanpad{3} = zeros([stack_size(1)+nanpix(1)*2 stack_size(2)+nanpix(2)*2 nanpix(3)]);
            tempImg = padnans(tempImg, nanpad);
            usave_out.multigrid_ops.nanpix = nanpix;
            clear nanpad;
        else
            udiff = diff(cat(3, stack_size(1:3), u_basis),[],3);
            if all(mod(udiff,2)==0) % this is a cheap and potentially wrong solution if nanpix isn't available
                nanpix = udiff/2;
                nanpad{1} = zeros([stack_size(1) nanpix(1) stack_size(3)]);
                nanpad{2} = zeros([nanpix(2) stack_size(2)+nanpix(1)*2 stack_size(3)]);
                nanpad{3} = zeros([stack_size(1)+nanpix(1)*2 stack_size(2)+nanpix(2)*2 nanpix(3)]);
                tempImg = padnans(tempImg, nanpad);
                usave_out.multigrid_ops.nanpix = nanpix;
                clear nanpad;
            else
                error; % not actually an error, but I've not yet written in the solution (JPM)
            end
        end
    else
        correct = false;
    end
    
    thisReg = register_multigrid_warp(tempImg,thisu,usave_out.rmg);
    if correct == true
        expsz = size(thisReg);
        thisReg = thisReg.*usave_out.rmg.image_grid(1).mask;
        thisReg = thisReg(1+nanpix(1):expsz(1)-nanpix(1), 1+nanpix(2):expsz(2)-nanpix(2), 1+nanpix(3):expsz(3)-nanpix(3));
    else
        thisReg = thisReg.*usave_out.rmg.image_grid(1).mask;
    end
    if strcmp('DV8285_BV', header.camera) 
      thisReg = cat(1,thisReg,padmat);
    end
    count = fwrite(fid,uint16(thisReg),'uint16');
    if (count < prod(stack_size(1:end-1)))
        error('Wrote fewer pixels than expected, quitting');
    end
    fprintf('%d..',idx);
    if mod(idx,10)==0
        fprintf('\n');
    end
end