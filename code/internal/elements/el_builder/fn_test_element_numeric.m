function [el_K, el_M, el_C, time_taken] = fn_test_element_numeric(el_fname, nds, no_trials)
%Function to test element function file to generate element matrices for
%specified number of elements. only one set of nodal coordinates is
%required; if no_trials>1, nodal coordinates are generated automatically by
%perturbing those specified. 
rng(1);
els = 1:size(nds, 1);
if no_trials > 1
    nds = repmat(nds, [no_trials, 1]);
    nds = nds + randn(size(nds)) * 0.1;
    els = (1:no_trials)' + els;
end

%In case filename is passed, extract function name
[~ ,fn_name] = fileparts(el_fname);

fn_el_type = str2func(fn_name);
% fn_el_mats_ref = str2func(['fn_el_', ref_element]);

%Call with empty arguments to get loc_nds and loc_dfs only
[~, ~, ~, loc_nds, loc_dfs] = fn_el_type([], [], [], []);
if any(loc_dfs == 4)
    solid_or_fluid = 'fluid';
    D = rand(1);
else
    solid_or_fluid = 'solid';
    D = rand(6);
    D = D + D';
end
rho = rand(1);
%--------------------------------------------------------------------------
%Do the test
tic;
[el_K, el_C, el_M, loc_nd, loc_df] = fn_el_type(nds, els, D, rho);
time_taken = double(toc);
% 
% if no_trials == 1
%     if nargout == 0
%         fprintf('test_el_K = \n');
%         disp(squeeze(test_el_K));
%         fprintf('test_el_M = \n');
%         disp(squeeze(test_el_M));
%     end
%     el_K = squeeze(test_el_K);
%     el_M = squeeze(test_el_M);
% else
%     el_K = [];
%     el_M = [];
% end
% 
% if exist(func2str(fn_el_mats_ref), 'file')
%     %Existing function
%     tic;
%     [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats_ref(nds, els, D, rho);
%     t_ref = double(toc);
%     if no_trials == 1
%         if nargout == 0
%             fprintf('el_K = \n');
%             disp(squeeze(el_K));
%             fprintf('el_M = \n');
%             disp(squeeze(el_M));
%         end
%         el_K_ref = squeeze(el_K);
%         el_M_ref = squeeze(el_M);
% 
%     else
%         el_K_ref = [];
%         el_M_ref = [];
%     end
%     %Comparison
%     if nargout == 0
%         fprintf(['\nCOMPARISON OF OUTPUTS FROM  ',func2str(fn_el_type),' AND ',func2str(fn_el_mats_ref),':\n'])
%         fprintf('  Fractional RMS error for K: %e\n', fn_compare_matrices(test_el_K, el_K));
%         fprintf('  Fractional RMS error for M: %e\n', fn_compare_matrices(test_el_M, el_M));
%     end
% else
%     t_ref = -1;
%     el_K_ref = [];
% end
% 
% if nargout == 0
%     fprintf('\nSPEED TEST FOR %i ELEMENTS\n', no_trials);
%     fprintf(['  ', func2str(fn_el_type), ' took %.2fs\n'], t_test);
%     if t_ref >= 0
%         fprintf(['  ', func2str(fn_el_mats_ref), ' took %.2fs\n'], t_ref);
%     end
% end

end