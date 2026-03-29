function fn_el_test(test_element, ref_element, parent_nds, trial_pts)

fn_el_mats_test = str2func(['fn_el_', test_element]);
fn_el_mats_ref = str2func(['fn_el_', ref_element]);

[~, ~, ~, ~, loc_dfs] = fn_el_mats_test([], [], [], []);
if any(loc_dfs == 4)
    solid_or_fluid = 'fluid';
else
    solid_or_fluid = 'solid';
end

nds = repmat(parent_nds, [trial_pts, 1]);
els = reshape(1:size(nds, 1), size(parent_nds, 1), [])';

if trial_pts > 1
    nds = nds + randn(size(nds)) * 0.1;
end

if trial_pts > 1
    rho = rand(1);
    switch solid_or_fluid
        case 'solid'
            D = rand(6);
            D = D + D';
            D = eye(6);
        case 'fluid'
            if trial_pts > 1
                D = rand(1);
            end
    end
else
    rho = 1;
    switch solid_or_fluid
        case 'solid'
            D = eye(6);
        case 'fluid'
            D = 1;
    end
end

%--------------------------------------------------------------------------
%Do the test
tic;
[test_el_K, test_el_C, test_el_M, test_loc_nd, test_loc_df] = fn_el_mats_test(nds, els, D, rho);
t_test = double(toc);

if trial_pts == 1
    fprintf('test_el_K = \n')
    disp(squeeze(test_el_K));
    fprintf('test_el_M = \n')
    disp(squeeze(test_el_M));
end

if exist(func2str(fn_el_mats_ref), 'file')
    %Existing function
    tic;
    [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats_ref(nds, els, D, rho);
    t_ref = double(toc);
    if trial_pts == 1
        fprintf('el_K = \n')
        disp(squeeze(el_K));
        fprintf('el_M = \n')
        disp(squeeze(el_M));
    end
    %Comparison
    fprintf(['\nCOMPARISON OF OUTPUTS FROM  ',func2str(fn_el_mats_test),' AND ',func2str(fn_el_mats_ref),':\n'])
    fprintf('  Fractional RMS error for K: %e\n', fn_compare_matrices(test_el_K, el_K));
    fprintf('  Fractional RMS error for M: %e\n', fn_compare_matrices(test_el_M, el_M));
else
    t_ref = -1;
end

fprintf('\nSPEED TEST FOR %i ELEMENTS\n', trial_pts);
fprintf(['  ', func2str(fn_el_mats_test), ' took %.2fs\n'], t_test);
if t_ref >= 0
    fprintf(['  ', func2str(fn_el_mats_ref), ' took %.2fs\n'], t_ref);
end

end