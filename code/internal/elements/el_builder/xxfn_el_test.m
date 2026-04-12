function [el_K_test, el_M_test, el_K_ref, el_M_ref] = fn_el_test(test_element, ref_element, solid_or_fluid, test_nds, test_D, test_rho, trial_pts)

els = 1:size(test_nds, 1);
if trial_pts > 1
    test_nds = repmat(test_nds, [trial_pts, 1]);
    test_nds = test_nds + randn(size(test_nds)) * 0.1;
    els = (1:trial_pts)' + els;
end

fn_el_mats_test = str2func(['fn_el_', test_element]);
fn_el_mats_ref = str2func(['fn_el_', ref_element]);

[~, ~, ~, ~, loc_dfs] = fn_el_mats_test([], [], [], []);
if any(loc_dfs == 4)
    solid_or_fluid = 'fluid';
else
    solid_or_fluid = 'solid';
end

% test_nds = repmat(test_nds, [trial_pts, 1]);

%--------------------------------------------------------------------------
%Do the test
tic;
[test_el_K, test_el_C, test_el_M, test_loc_nd, test_loc_df] = fn_el_mats_test(test_nds, els, test_D, test_rho);
t_test = double(toc);

if trial_pts == 1
    if nargout == 0
        fprintf('test_el_K = \n');
        disp(squeeze(test_el_K));
        fprintf('test_el_M = \n');
        disp(squeeze(test_el_M));
    end
    el_K_test = squeeze(test_el_K);
    el_M_test = squeeze(test_el_M);
else
    el_K_test = [];
    el_M_test = [];
end

if exist(func2str(fn_el_mats_ref), 'file')
    %Existing function
    tic;
    [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats_ref(test_nds, els, test_D, test_rho);
    t_ref = double(toc);
    if trial_pts == 1
        if nargout == 0
            fprintf('el_K = \n');
            disp(squeeze(el_K));
            fprintf('el_M = \n');
            disp(squeeze(el_M));
        end
        el_K_ref = squeeze(el_K);
        el_M_ref = squeeze(el_M);

    else
        el_K_ref = [];
        el_M_ref = [];
    end
    %Comparison
    if nargout == 0
        fprintf(['\nCOMPARISON OF OUTPUTS FROM  ',func2str(fn_el_mats_test),' AND ',func2str(fn_el_mats_ref),':\n'])
        fprintf('  Fractional RMS error for K: %e\n', fn_compare_matrices(test_el_K, el_K));
        fprintf('  Fractional RMS error for M: %e\n', fn_compare_matrices(test_el_M, el_M));
    end
else
    t_ref = -1;
    el_K_ref = [];
end

if nargout == 0
    fprintf('\nSPEED TEST FOR %i ELEMENTS\n', trial_pts);
    fprintf(['  ', func2str(fn_el_mats_test), ' took %.2fs\n'], t_test);
    if t_ref >= 0
        fprintf(['  ', func2str(fn_el_mats_ref), ' took %.2fs\n'], t_ref);
    end
end

end