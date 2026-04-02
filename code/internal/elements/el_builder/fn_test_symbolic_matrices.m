function [K, M] = fn_test_symbolic_matrices(sym_mats, test_nds, test_D)
%Purpose of this is to allow tests of the symbolic matrices before they are
%converted into m-files, so that errors in the symbolic formulations can be
%identified first

no_spatial_dims = 3;
K = zeros(numel(sym_mats.loc_nd));
no_dfs = numel(unique(sym_mats.loc_df));
no_dims_of_el = size(test_nds, 2);
for g = 1:numel(sym_mats.gauss_wts)
    detJ = subs(sym_mats.detJ(g), sym_mats.nds, test_nds);
    if isfield(sym_mats, 'K')
        Kg = sym_mats.K(:,:,g);
        Kg = subs(Kg, sym_mats.nds, test_nds);
        Kg = subs(Kg, 'detJ', detJ);
    else
        if isfield(sym_mats, 'B')
            B = sym_mats.B(:,:,g);
            B = subs(B, sym_mats.nds, test_nds);
            B = subs(B, 'detJ', detJ);
        else
            N_diff = sym_mats.N_diff(:, :, g);
            E = zeros(no_dfs * no_dims_of_el, numel(sym_mats.loc_nd));
            for n = 1:size(N_diff, 2) %outer loop across columns in steps of no_dfs
                for j = 1:no_dfs %outer loop down rows (first strain index - disp comp)
                    for k = 1:size(N_diff, 1) %inner loop down rows (second strain index - deriv direction)
                        E((j - 1) * no_dims_of_el + k, (n-1) * no_dfs + j) = N_diff(k, n);
                    end
                end
            end
            invJ = sym_mats.invJ(:, :, g);
            invJ = subs(invJ, sym_mats.nds, test_nds);
            invJ = subs(invJ, 'detJ', detJ);
            invJpadded = zeros(no_spatial_dims, no_dims_of_el);
            invJpadded(1:no_dims_of_el, :) = invJ;
            invJstar = kron(eye(no_dfs), invJpadded);
            B = sym_mats.L * invJstar * E;
        end
        Kg = B' * test_D *  B * detJ;
    end
    K = K + Kg * sym_mats.gauss_wts(g);
end
M = [];
end