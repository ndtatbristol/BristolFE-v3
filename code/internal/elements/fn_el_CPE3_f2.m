function [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_CPE3_f2(nds, els, D, rho, varargin)
%SUMMARY
%	This function was created automatically by fn_create_element_matrix_file3
%	and contains code to return the stiffness and mass matrices
%	for multiple elements of the same material and type given by the latter
%	part of the filename, fn_el_CPE3_f2.
%INPUTS
%	nds - n_nds x n_dims matrix of nodal coordinates
%	els - n_els x n_nds_per_el matrix of node indices for each elements
%	D - ns x ns material stiffness matrix
%	rho - material density
%	[dofs_to_use = [] - optional vector listing the DoFs to use, e.g. [1, 2]. Use [] for all]
%OUTPUTS
%	el_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices
%AUTHOR
%	Paul Wilcox (12-Apr-2026 09:06:00)

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
gauss_wts = zeros(1, 1);
gauss_wts(1) = 5.000000000000000000e-01;

%Zero the outputs
el_K = zeros(9, 9, no_els);
el_M_tmp = zeros(9, 9, no_els);
el_C = zeros(9, 9, no_els);

detJ = zeros(1, 1, no_els);
N = zeros(3, 9, no_els);
J = zeros(2, 2, no_els);
B = zeros(6, 9, no_els);
%Loop over Gauss points
for g = 1:no_gauss_pts

	switch g
		case 1
			%Terms of Jacobian
			J(1, 1, :) = nds_2_1 - nds_1_1;
			J(1, 2, :) = nds_2_2 - nds_1_2;
			J(2, 1, :) = nds_3_1 - nds_1_1;
			J(2, 2, :) = nds_3_2 - nds_1_2;

			%Determinant of Jacobian
			detJ = J(1, 1, :) .* J(2, 2, :) - J(1, 2, :) .* J(2, 1, :);

			%Terms of B matrix
			B(1, 1, :) = J(1, 2, :) ./ detJ - J(2, 2, :) ./ detJ;
			B(1, 4, :) = J(2, 2, :) ./ detJ;
			B(1, 7, :) = -J(1, 2, :) ./ detJ;
			B(2, 2, :) = J(2, 1, :) ./ detJ - J(1, 1, :) ./ detJ;
			B(2, 5, :) = -J(2, 1, :) ./ detJ;
			B(2, 8, :) = J(1, 1, :) ./ detJ;
			B(4, 3, :) = J(2, 1, :) ./ detJ - J(1, 1, :) ./ detJ;
			B(4, 6, :) = -J(2, 1, :) ./ detJ;
			B(4, 9, :) = J(1, 1, :) ./ detJ;
			B(5, 3, :) = J(1, 2, :) ./ detJ - J(2, 2, :) ./ detJ;
			B(5, 6, :) = J(2, 2, :) ./ detJ;
			B(5, 9, :) = -J(1, 2, :) ./ detJ;
			B(6, 1, :) = J(2, 1, :) ./ detJ - J(1, 1, :) ./ detJ;
			B(6, 2, :) = J(1, 2, :) ./ detJ - J(2, 2, :) ./ detJ;
			B(6, 4, :) = -J(2, 1, :) ./ detJ;
			B(6, 5, :) = J(2, 2, :) ./ detJ;
			B(6, 7, :) = J(1, 1, :) ./ detJ;
			B(6, 8, :) = -J(1, 2, :) ./ detJ;

			%Terms of N matrix
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

	%Evaluate K = B'DB|J| and accumulate over Gauss points
	el_K = el_K + pagemtimes(pagemtimes(B, 'transpose', pagemtimes(D, B), 'none'), detJ) * gauss_wts(g);

	%Evaluate rho * N'N|J|
	el_M_tmp = el_M_tmp + rho * pagemtimes(pagemtimes(N, 'transpose', N, 'none'), detJ) * gauss_wts(g);

end

%Diagonalise M
el_M = zeros(9, 9, no_els);
for i = 1:9
	el_M(i, i, :) = sum(el_M_tmp(:, i, :), 1);
end

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
