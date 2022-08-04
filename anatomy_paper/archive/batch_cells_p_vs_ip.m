%%% run chi square test for ever row in a table - must already have patchdata, and cells image in rows of filetable
% f is filetable
%%%% updated 18/9/19

use_cell_dense_roi = 0; % use the cell-dense roi rather than the full roi from patchdata
do_permtest = 1;
    permutation_pars = struct;
    permutation_pars.npermutations = 10; % number of permutations; ?10,000 recommended to asymptote
    permutation_pars.custom_interpatches = 1; % use something other than all non-patch pixels as interpatch pixels
    permutation_pars.interpatchQuantiles = [1 ]; % quantiles selected to be counted as interpatches, if custom_interpatches==true; lower numbers = darker pixels
do_chitest = 1;

for i = 1:height(f)
    if do_permtest
        if ~ismember('permtest',f.Properties.VariableNames)
            f.permtest = cell(height(f),1);
        end
        if use_cell_dense_roi
            patchimage_cell_dense_roi = patchdata.patchimage & loadbw(f.cell_dense_roifile{i});
            f.permtest{i} = permutation_test_patches_cells(patchimage_cell_dense_roi,f.cellimage{i},f.cell_dense_roifile{i},f.zoom(i),patchdata.scope,permutation_pars);
        else
            f.permtest{i} = permutation_test_patches_cells(patchdata.patchimage,f.cellimage{i},patchdata.roi,f.zoom(i),patchdata.scope,permutation_pars);
        end
        f.patchInterpatchRatio(i) = f.permtest{i}.patchInterpatchRatio;
        f.pval_permutation_test(i) = f.permtest{i}.pval_permutation_test;
        f.ipatch_cellsPerSqMM(i) = f.permtest{i}.interpatch_cellsPerSqMM;
    end
    if do_chitest
        if use_cell_dense_roi
            [f.chi_p(i), f.quantcells{i}, f.cellcenters{i}] = chi_square_quantiles_cells(patchdata, f.cellimage{i}, f.cell_dense_roifile{i});
        else
            [f.chi_p(i), f.quantcells{i}, f.cellcenters{i}] = chi_square_quantiles_cells(patchdata, f.cellimage{i});
        end
    end
end