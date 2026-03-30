function [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_CPE3_new(nds, els, D, rho, varargin)
%SUMMARY
%	This function was created automatically by fn_create_element_matrix_file3
%	and contains code to return the stiffness and mass matrices
%	for multiple elements of the same material and type given by the latter
%	part of the filename, fn_el_CPE3_new.
%INPUTS
%	nds - n_nds x n_dims matrix of nodal coordinates
%	els - n_els x n_nds_per_el matrix of node indices for each elements
%	D - ns x ns material stiffness matrix
%	rho - material density
%	[dofs_to_use = [] - optional vector listing the DoFs to use, e.g. [1, 2]. Use [] for all]
%OUTPUTS
%	el_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices
%AUTHOR
%	Paul Wilcox (30-Mar-2026 13:48:48)

%Define sqrt(3)
rt3 = sqrt(3);
%Deal with optional argument about which DOFs to use
if isempty(varargin)
	dofs_to_use = [];
else
	dofs_to_use = varargin{1};
end

%Record the local node numbers of the element stiffness matrices
loc_nd = [1  1  1  2  2  2  3  3  3];

%Record the local DOFs of the element stiffness matrices
loc_df = [1  2  3  1  2  3  1  2  3];

%Get the DOFs if not specified
if isempty(dofs_to_use)
	dofs_to_use = unique(loc_df);
end

%If any inputs blank, return at this point with just the loc_nd and loc_df
if isempty(nds) || isempty(els) || isempty(D) || isempty(rho)
	el_K = [];
	el_M = [];
	el_C = [];
	[loc_nd, loc_df] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use);
	return
end


%Some constants
no_gauss_pts = 1;
no_els = size(els, 1);

%Matrices of nodal coordinates
nds_1_1 = nds(els(:, 1), 1);
nds_1_2 = nds(els(:, 1), 2);
nds_2_1 = nds(els(:, 2), 1);
nds_2_2 = nds(els(:, 2), 2);
nds_3_1 = nds(els(:, 3), 1);
nds_3_2 = nds(els(:, 3), 2);

%Vector of Gauss weights
W = zeros(1, 1);
W(1) = 5.000000000000000000e-01;

%Zero the outputs
el_K = zeros(9, 9, no_els);
el_M_tmp = zeros(9, 9, no_els);
el_C = zeros(9, 9, no_els);

detJ = zeros(no_els, 1);
B = zeros(6, 9, no_els);
N = zeros(3, 9, no_els);
%Loop over Gauss points
for i = 1:no_gauss_pts

    %Jacobians, N- and B-matrices at each Gauss point

    switch i
        case 1
            detJ = nds_1_1 .* nds_2_2 - nds_1_2 .* nds_2_1 - nds_1_1 .* nds_3_2 + nds_1_2 .* nds_3_1 + nds_2_1 .* nds_3_2 - nds_2_2 .* nds_3_1;

            B(1, 1, :) = (nds_2_2 - nds_3_2) ./ detJ;
            B(1, 4, :) = -(nds_1_2 - nds_3_2) ./ detJ;
            B(1, 7, :) = (nds_1_2 - nds_2_2) ./ detJ;
            B(2, 2, :) = -(nds_2_1 - nds_3_1) ./ detJ;
            B(2, 5, :) = (nds_1_1 - nds_3_1) ./ detJ;
            B(2, 8, :) = -(nds_1_1 - nds_2_1) ./ detJ;
            B(4, 3, :) = -(nds_2_1 - nds_3_1) ./ detJ;
            B(4, 6, :) = (nds_1_1 - nds_3_1) ./ detJ;
            B(4, 9, :) = -(nds_1_1 - nds_2_1) ./ detJ;
            B(5, 3, :) = (nds_2_2 - nds_3_2) ./ detJ;
            B(5, 6, :) = -(nds_1_2 - nds_3_2) ./ detJ;
            B(5, 9, :) = (nds_1_2 - nds_2_2) ./ detJ;
            B(6, 1, :) = -(nds_2_1 - nds_3_1) ./ detJ;
            B(6, 2, :) = (nds_2_2 - nds_3_2) ./ detJ;
            B(6, 4, :) = (nds_1_1 - nds_3_1) ./ detJ;
            B(6, 5, :) = -(nds_1_2 - nds_3_2) ./ detJ;
            B(6, 7, :) = -(nds_1_1 - nds_2_1) ./ detJ;
            B(6, 8, :) = (nds_1_2 - nds_2_2) ./ detJ;

            N(1, 1, :) = 1 ./ 3;
            N(1, 4, :) = 1 ./ 3;
            N(1, 7, :) = 1 ./ 3;
            N(2, 2, :) = 1 ./ 3;
            N(2, 5, :) = 1 ./ 3;
            N(2, 8, :) = 1 ./ 3;
            N(3, 3, :) = 1 ./ 3;
            N(3, 6, :) = 1 ./ 3;
            N(3, 9, :) = 1 ./ 3;

    end

    %Evaluate B'DB|J|
    el_K = el_K + pagemtimes(pagemtimes(B, 'transpose', pagemtimes(D, B), 'none'), permute(detJ, [2, 3, 1])) * W(i);


    %Evaluate N'rhoN|J|
    el_M_tmp = el_M_tmp + pagemtimes(pagemtimes(N, 'transpose', rho * N, 'none'), permute(detJ, [2, 3, 1])) * W(i);

end

%Diagonalise M
el_M = zeros(9, 9, no_els);
el_M(1, 1, :) = sum(el_M_tmp(:, 1, :), 1);
el_M(2, 2, :) = sum(el_M_tmp(:, 2, :), 1);
el_M(3, 3, :) = sum(el_M_tmp(:, 3, :), 1);
el_M(4, 4, :) = sum(el_M_tmp(:, 4, :), 1);
el_M(5, 5, :) = sum(el_M_tmp(:, 5, :), 1);
el_M(6, 6, :) = sum(el_M_tmp(:, 6, :), 1);
el_M(7, 7, :) = sum(el_M_tmp(:, 7, :), 1);
el_M(8, 8, :) = sum(el_M_tmp(:, 8, :), 1);
el_M(9, 9, :) = sum(el_M_tmp(:, 9, :), 1);

%Remove unwanted DOFs from element matrices
j = ismember(loc_df, dofs_to_use);
el_K = el_K(j, j, :);
el_M = el_M(j, j, :);
el_C = el_C(j, j, :);
loc_nd = loc_nd(j);
loc_df = loc_df(j);

%Change dimension order of element matrices
el_K = permute(el_K, [3, 1, 2]);
el_M = permute(el_M, [3, 1, 2]);
el_C = permute(el_C, [3, 1, 2]);

end
