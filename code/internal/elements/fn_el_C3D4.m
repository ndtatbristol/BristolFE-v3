function [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_C3D4(nds, els, D, rho, varargin)
%SUMMARY
%	This function was created automatically by fn_create_element_matrix_file3
%	and contains code to return the stiffness and mass matrices
%	for multiple elements of the same material and type given by the latter
%	part of the filename, fn_el_C3D4.
%INPUTS
%	nds - n_nds x n_dims matrix of nodal coordinates
%	els - n_els x n_nds_per_el matrix of node indices for each elements
%	D - ns x ns material stiffness matrix
%	rho - material density
%	[dofs_to_use = [] - optional vector listing the DoFs to use, e.g. [1, 2]. Use [] for all]
%OUTPUTS
%	el_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices
%AUTHOR
%	Paul Wilcox (09-Jun-2026 10:54:06)

%--------------------------------------------------------------------------

%Deal with optional argument about which DOFs to use
if numel(varargin) < 1
	dofs_to_use = [];
else
	dofs_to_use = varargin{1};
end
if numel(varargin) < 2
	use_gm_builder_v5_dim_order = 1;
else
	use_gm_builder_v5_dim_order = varargin{2};
end

%Record the local node numbers of the element stiffness matrices
loc_nd = [1  1  1  2  2  2  3  3  3  4  4  4];

%Record the local DOFs of the element stiffness matrices
loc_df = [1  2  3  1  2  3  1  2  3  1  2  3];

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

%Constants
no_gauss_pts = 1;
no_els = size(els, 1);
root3 = sqrt(3);

%Matrices of nodal coordinates
nds_1_1 = nds(els(:, 1), 1);
nds_1_2 = nds(els(:, 1), 2);
nds_1_3 = nds(els(:, 1), 3);
nds_2_1 = nds(els(:, 2), 1);
nds_2_2 = nds(els(:, 2), 2);
nds_2_3 = nds(els(:, 2), 3);
nds_3_1 = nds(els(:, 3), 1);
nds_3_2 = nds(els(:, 3), 2);
nds_3_3 = nds(els(:, 3), 3);
nds_4_1 = nds(els(:, 4), 1);
nds_4_2 = nds(els(:, 4), 2);
nds_4_3 = nds(els(:, 4), 3);

%Vector of Gauss weights
gauss_wts = zeros(1, 1);
gauss_wts(1) = 1.666666666666666574e-01;

%Zero the outputs
el_K = zeros(12, 12, no_els);
el_C = zeros(12, 12, no_els);
el_M_tmp = zeros(12, 12, no_els);
detJ = zeros(1, 1, no_els);
N = zeros(3, 12, no_els);
J = zeros(3, 3, no_els);

B2 = zeros(9, 9, no_els);
B3 = zeros(9, 12);
%Factors of B matrix are B1, B2, and B3. Only B2 is a function of the specific
 %element. B1 is also independent of Gauss point and is defined first.
B1 = [
	1, 0, 0, 0, 0, 0, 0, 0, 0
	0, 0, 0, 0, 1, 0, 0, 0, 0
	0, 0, 0, 0, 0, 0, 0, 0, 1
	0, 0, 0, 0, 0, 1, 0, 1, 0
	0, 0, 1, 0, 0, 0, 1, 0, 0
	0, 1, 0, 1, 0, 0, 0, 0, 0
];
%Loop over Gauss points
for g = 1:no_gauss_pts

	switch g
	%Define matrices that depend on Gauss point
		case 1
			%Terms of Jacobian
			J(1, 1, :) = nds_2_1 - nds_1_1;
			J(1, 2, :) = nds_2_2 - nds_1_2;
			J(1, 3, :) = nds_2_3 - nds_1_3;
			J(2, 1, :) = nds_3_1 - nds_1_1;
			J(2, 2, :) = nds_3_2 - nds_1_2;
			J(2, 3, :) = nds_3_3 - nds_1_3;
			J(3, 1, :) = nds_4_1 - nds_1_1;
			J(3, 2, :) = nds_4_2 - nds_1_2;
			J(3, 3, :) = nds_4_3 - nds_1_3;

			%Terms of B3 matrix
			B3(1, 1, :) = -1;
			B3(1, 4, :) = 1;
			B3(2, 1, :) = -1;
			B3(2, 7, :) = 1;
			B3(3, 1, :) = -1;
			B3(3, 10, :) = 1;
			B3(4, 2, :) = -1;
			B3(4, 5, :) = 1;
			B3(5, 2, :) = -1;
			B3(5, 8, :) = 1;
			B3(6, 2, :) = -1;
			B3(6, 11, :) = 1;
			B3(7, 3, :) = -1;
			B3(7, 6, :) = 1;
			B3(8, 3, :) = -1;
			B3(8, 9, :) = 1;
			B3(9, 3, :) = -1;
			B3(9, 12, :) = 1;

			%Terms of N matrix
			N(1, 1, :) = 1 ./ 4;
			N(1, 4, :) = 1 ./ 4;
			N(1, 7, :) = 1 ./ 4;
			N(1, 10, :) = 1 ./ 4;
			N(2, 2, :) = 1 ./ 4;
			N(2, 5, :) = 1 ./ 4;
			N(2, 8, :) = 1 ./ 4;
			N(2, 11, :) = 1 ./ 4;
			N(3, 3, :) = 1 ./ 4;
			N(3, 6, :) = 1 ./ 4;
			N(3, 9, :) = 1 ./ 4;
			N(3, 12, :) = 1 ./ 4;

	end

	%Determinant of Jacobian
	detJ = J(1, 1, :) .* J(2, 2, :) .* J(3, 3, :) - J(1, 1, :) .* J(2, 3, :) .* J(3, 2, :) - J(1, 2, :) .* J(2, 1, :) .* J(3, 3, :) + J(1, 2, :) .* J(2, 3, :) .* J(3, 1, :) + J(1, 3, :) .* J(2, 1, :) .* J(3, 2, :) - J(1, 3, :) .* J(2, 2, :) .* J(3, 1, :);

	%Terms of B2 matrix
	B2(1, 1, :) = (J(2, 2, :) .* J(3, 3, :) - J(2, 3, :) .* J(3, 2, :)) ./ detJ;
	B2(1, 2, :) = -(J(1, 2, :) .* J(3, 3, :) - J(1, 3, :) .* J(3, 2, :)) ./ detJ;
	B2(1, 3, :) = (J(1, 2, :) .* J(2, 3, :) - J(1, 3, :) .* J(2, 2, :)) ./ detJ;
	B2(2, 1, :) = -(J(2, 1, :) .* J(3, 3, :) - J(2, 3, :) .* J(3, 1, :)) ./ detJ;
	B2(2, 2, :) = (J(1, 1, :) .* J(3, 3, :) - J(1, 3, :) .* J(3, 1, :)) ./ detJ;
	B2(2, 3, :) = -(J(1, 1, :) .* J(2, 3, :) - J(1, 3, :) .* J(2, 1, :)) ./ detJ;
	B2(3, 1, :) = (J(2, 1, :) .* J(3, 2, :) - J(2, 2, :) .* J(3, 1, :)) ./ detJ;
	B2(3, 2, :) = -(J(1, 1, :) .* J(3, 2, :) - J(1, 2, :) .* J(3, 1, :)) ./ detJ;
	B2(3, 3, :) = (J(1, 1, :) .* J(2, 2, :) - J(1, 2, :) .* J(2, 1, :)) ./ detJ;
	B2(4, 4, :) = (J(2, 2, :) .* J(3, 3, :) - J(2, 3, :) .* J(3, 2, :)) ./ detJ;
	B2(4, 5, :) = -(J(1, 2, :) .* J(3, 3, :) - J(1, 3, :) .* J(3, 2, :)) ./ detJ;
	B2(4, 6, :) = (J(1, 2, :) .* J(2, 3, :) - J(1, 3, :) .* J(2, 2, :)) ./ detJ;
	B2(5, 4, :) = -(J(2, 1, :) .* J(3, 3, :) - J(2, 3, :) .* J(3, 1, :)) ./ detJ;
	B2(5, 5, :) = (J(1, 1, :) .* J(3, 3, :) - J(1, 3, :) .* J(3, 1, :)) ./ detJ;
	B2(5, 6, :) = -(J(1, 1, :) .* J(2, 3, :) - J(1, 3, :) .* J(2, 1, :)) ./ detJ;
	B2(6, 4, :) = (J(2, 1, :) .* J(3, 2, :) - J(2, 2, :) .* J(3, 1, :)) ./ detJ;
	B2(6, 5, :) = -(J(1, 1, :) .* J(3, 2, :) - J(1, 2, :) .* J(3, 1, :)) ./ detJ;
	B2(6, 6, :) = (J(1, 1, :) .* J(2, 2, :) - J(1, 2, :) .* J(2, 1, :)) ./ detJ;
	B2(7, 7, :) = (J(2, 2, :) .* J(3, 3, :) - J(2, 3, :) .* J(3, 2, :)) ./ detJ;
	B2(7, 8, :) = -(J(1, 2, :) .* J(3, 3, :) - J(1, 3, :) .* J(3, 2, :)) ./ detJ;
	B2(7, 9, :) = (J(1, 2, :) .* J(2, 3, :) - J(1, 3, :) .* J(2, 2, :)) ./ detJ;
	B2(8, 7, :) = -(J(2, 1, :) .* J(3, 3, :) - J(2, 3, :) .* J(3, 1, :)) ./ detJ;
	B2(8, 8, :) = (J(1, 1, :) .* J(3, 3, :) - J(1, 3, :) .* J(3, 1, :)) ./ detJ;
	B2(8, 9, :) = -(J(1, 1, :) .* J(2, 3, :) - J(1, 3, :) .* J(2, 1, :)) ./ detJ;
	B2(9, 7, :) = (J(2, 1, :) .* J(3, 2, :) - J(2, 2, :) .* J(3, 1, :)) ./ detJ;
	B2(9, 8, :) = -(J(1, 1, :) .* J(3, 2, :) - J(1, 2, :) .* J(3, 1, :)) ./ detJ;
	B2(9, 9, :) = (J(1, 1, :) .* J(2, 2, :) - J(1, 2, :) .* J(2, 1, :)) ./ detJ;

	%Calculate B matrix
	B = pagemtimes(B1, pagemtimes(B2, B3));

	%Evaluate K = B'DB|J| and accumulate over Gauss points
	el_K = el_K + pagemtimes(pagemtimes(B, 'transpose', pagemtimes(D, B), 'none'), detJ) * gauss_wts(g);

	%Evaluate rho * N'N|J|
	el_M_tmp = el_M_tmp + rho * pagemtimes(pagemtimes(N, 'transpose', N, 'none'), detJ) * gauss_wts(g);

end

%Diagonalise M
el_M = zeros(12, 12, no_els);
for i = 1:12
	el_M(i, i, :) = sum(el_M_tmp(:, i, :), 1);
end

%Remove unwanted DOFs from element matrices
j = ismember(loc_df, dofs_to_use);
el_K = el_K(j, j, :);
el_M = el_M(j, j, :);
el_C = el_C(j, j, :);
loc_nd = loc_nd(j);
loc_df = loc_df(j);

%Change dimension order of element matrices for v5 global matrix builder
if use_gm_builder_v5_dim_order
	el_K = permute(el_K, [3, 1, 2]);
	el_M = permute(el_M, [3, 1, 2]);
	el_C = permute(el_C, [3, 1, 2]);
end

end
