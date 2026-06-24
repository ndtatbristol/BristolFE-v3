function [K, M] = fn_test_symbolic_matrices(sym_mats, test_nds, test_D, test_rho)
%Purpose of this is to allow tests of the symbolic matrices before they are
%converted into m-files, so that errors in the symbolic formulations can be
%identified first
no_gauss_pts = numel(sym_mats.gauss_wts);
no_dfs = numel(unique(sym_mats.loc_df));
no_dims_of_el = size(test_nds, 2);

%insert - get form of symbolic K for testing
% B = sym_mats.B1 * sym_mats.B2 * sym_mats.B3;
% S = sym('J_%d_%d', size(sym_mats.J));
% B = subs(B, S, sym_mats.J);
% K = B' * sym('rho', 'real') * B * sym('detJ', 'real') * sym_mats.gauss_wts(1);
% 
% disp(char(simplify(K(1,1), 10)));

B1 = sym_mats.B1;
if isscalar(test_D)
    scaling = double(subs(subs(sym_mats.scaling, sym_mats.rho, test_rho), sym_mats.D, test_D));
else
    scaling = 1;
end

K = zeros(numel(sym_mats.loc_nd));
M = zeros(numel(sym_mats.loc_nd));
for g = 1:no_gauss_pts
    J = subs(sym_mats.J(:,:,g), sym_mats.nds_sym, test_nds);
    detJ = double(subs(sym_mats.detJ, sym_mats.J_sym, J));
    B2 = double(subs(subs(sym_mats.B2, sym_mats.J_sym, J), sym_mats.detJ_sym, detJ));
    B3 = double(sym_mats.B3(:, :, g));
    N = double(sym_mats.N(:, :, g));
    B = B1 * B2 * B3;
    K = K + B' * test_D * B * detJ * sym_mats.gauss_wts(g);
    M = M + N' * test_rho * N * detJ * sym_mats.gauss_wts(g);
end
K = K * scaling;
M = M * scaling;

M = diag(sum(M));
    % if isfield(sym_mats, 'K')
    %     Kg = sym_mats.K(:,:,g);
    %     Kg = subs(Kg, sym_mats.nds, test_nds);
    %     Kg = subs(Kg, 'detJ', detJ);
    %     Kg = subs(Kg, sym_mats.D, test_D);
    % else
    %     if isfield(sym_mats, 'B')
    %         B = sym_mats.B(:,:,g);
    %         B = subs(B, sym_mats.nds, test_nds);
    %         B = subs(B, 'detJ', detJ);
    %     else
    %         N_diff = sym_mats.N_diff(:, :, g);
    %         E = zeros(no_dfs * no_dims_of_el, numel(sym_mats.loc_nd));
    %         for n = 1:size(N_diff, 2) %outer loop across columns in steps of no_dfs
    %             for j = 1:no_dfs %outer loop down rows (first strain index - disp comp)
    %                 for k = 1:size(N_diff, 1) %inner loop down rows (second strain index - deriv direction)
    %                     E((j - 1) * no_dims_of_el + k, (n-1) * no_dfs + j) = N_diff(k, n);
    %                 end
    %             end
    %         end
    %         invJ = sym_mats.invJ(:, :, g);
    %         invJ = subs(invJ, sym_mats.nds, test_nds);
    %         invJ = subs(invJ, 'detJ', detJ);
    %         invJpadded = zeros(no_spatial_dims, no_dims_of_el);
    %         invJpadded(1:no_dims_of_el, :) = invJ;
    %         invJstar = kron(eye(no_dfs), invJpadded);
    %         B = sym_mats.L * invJstar * E;
    %     end
    %     Kg = B' * test_D *  B * detJ;
    % end
%     K = K + double(Kg * sym_mats.gauss_wts(g));
% end
end