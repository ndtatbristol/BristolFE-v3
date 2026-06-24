function [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_CPE4(nds, els, D, rho, varargin)
%SUMMARY
%	This function was created automatically by fn_create_element_matrix_file3
%	and contains code to return the stiffness and mass matrices
%	for multiple elements of the same material and type given by the latter
%	part of the filename, fn_el_CPE4.
%INPUTS
%	nds - n_nds x n_dims matrix of nodal coordinates
%	els - n_els x n_nds_per_el matrix of node indices for each elements
%	D - ns x ns material stiffness matrix
%	rho - material density
%	[dofs_to_use = [] - optional vector listing the DoFs to use, e.g. [1, 2]. Use [] for all]
%OUTPUTS
%	el_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices
%AUTHOR
%	Paul Wilcox (05-Jun-2026 15:19:33)

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
no_gauss_pts = 4;
no_els = size(els, 1);
root3 = sqrt(3);

%Matrices of nodal coordinates
nds_1_1 = nds(els(:, 1), 1);
nds_1_2 = nds(els(:, 1), 2);
nds_2_1 = nds(els(:, 2), 1);
nds_2_2 = nds(els(:, 2), 2);
nds_3_1 = nds(els(:, 3), 1);
nds_3_2 = nds(els(:, 3), 2);
nds_4_1 = nds(els(:, 4), 1);
nds_4_2 = nds(els(:, 4), 2);

%Vector of Gauss weights
gauss_wts = zeros(4, 1);
gauss_wts(1) = 1.000000000000000000e+00;
gauss_wts(2) = 1.000000000000000000e+00;
gauss_wts(3) = 1.000000000000000000e+00;
gauss_wts(4) = 1.000000000000000000e+00;

%Zero the outputs
el_K = zeros(12, 12, no_els);
el_C = zeros(12, 12, no_els);
el_M_tmp = zeros(12, 12, no_els);
detJ = zeros(1, 1, no_els);
N = zeros(3, 12, no_els);
J = zeros(2, 2, no_els);

B2 = zeros(9, 6, no_els);
B3 = zeros(6, 12);
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
			J(1, 1, :) = (nds_2_1 .* (root3 ./ 3 - 1)) ./ 4 - (nds_1_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_3_1 .* (root3 ./ 3 - 1)) ./ 4 + (nds_4_1 .* (root3 ./ 3 + 1)) ./ 4;
			J(1, 2, :) = (nds_2_2 .* (root3 ./ 3 - 1)) ./ 4 - (nds_1_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_3_2 .* (root3 ./ 3 - 1)) ./ 4 + (nds_4_2 .* (root3 ./ 3 + 1)) ./ 4;
			J(2, 1, :) = (nds_2_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_1_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_3_1 .* (root3 ./ 3 - 1)) ./ 4 + (nds_4_1 .* (root3 ./ 3 - 1)) ./ 4;
			J(2, 2, :) = (nds_2_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_1_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_3_2 .* (root3 ./ 3 - 1)) ./ 4 + (nds_4_2 .* (root3 ./ 3 - 1)) ./ 4;

			%Terms of B3 matrix
			B3(1, 1, :) = - root3 ./ 12 - 1 ./ 4;
			B3(1, 4, :) = root3 ./ 12 - 1 ./ 4;
			B3(1, 7, :) = 1 ./ 4 - root3 ./ 12;
			B3(1, 10, :) = root3 ./ 12 + 1 ./ 4;
			B3(2, 1, :) = - root3 ./ 12 - 1 ./ 4;
			B3(2, 4, :) = root3 ./ 12 + 1 ./ 4;
			B3(2, 7, :) = 1 ./ 4 - root3 ./ 12;
			B3(2, 10, :) = root3 ./ 12 - 1 ./ 4;
			B3(3, 2, :) = - root3 ./ 12 - 1 ./ 4;
			B3(3, 5, :) = root3 ./ 12 - 1 ./ 4;
			B3(3, 8, :) = 1 ./ 4 - root3 ./ 12;
			B3(3, 11, :) = root3 ./ 12 + 1 ./ 4;
			B3(4, 2, :) = - root3 ./ 12 - 1 ./ 4;
			B3(4, 5, :) = root3 ./ 12 + 1 ./ 4;
			B3(4, 8, :) = 1 ./ 4 - root3 ./ 12;
			B3(4, 11, :) = root3 ./ 12 - 1 ./ 4;
			B3(5, 3, :) = - root3 ./ 12 - 1 ./ 4;
			B3(5, 6, :) = root3 ./ 12 - 1 ./ 4;
			B3(5, 9, :) = 1 ./ 4 - root3 ./ 12;
			B3(5, 12, :) = root3 ./ 12 + 1 ./ 4;
			B3(6, 3, :) = - root3 ./ 12 - 1 ./ 4;
			B3(6, 6, :) = root3 ./ 12 + 1 ./ 4;
			B3(6, 9, :) = 1 ./ 4 - root3 ./ 12;
			B3(6, 12, :) = root3 ./ 12 - 1 ./ 4;

			%Terms of N matrix
			N(1, 1, :) = root3 ./ 6 + 1 ./ 3;
			N(1, 4, :) = 1 ./ 6;
			N(1, 7, :) = 1 ./ 3 - root3 ./ 6;
			N(1, 10, :) = 1 ./ 6;
			N(2, 2, :) = root3 ./ 6 + 1 ./ 3;
			N(2, 5, :) = 1 ./ 6;
			N(2, 8, :) = 1 ./ 3 - root3 ./ 6;
			N(2, 11, :) = 1 ./ 6;
			N(3, 3, :) = root3 ./ 6 + 1 ./ 3;
			N(3, 6, :) = 1 ./ 6;
			N(3, 9, :) = 1 ./ 3 - root3 ./ 6;
			N(3, 12, :) = 1 ./ 6;

	%Define matrices that depend on Gauss point
		case 2
			%Terms of Jacobian
			J(1, 1, :) = (nds_1_1 .* (root3 ./ 3 - 1)) ./ 4 - (nds_2_1 .* (root3 ./ 3 + 1)) ./ 4 + (nds_3_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_4_1 .* (root3 ./ 3 - 1)) ./ 4;
			J(1, 2, :) = (nds_1_2 .* (root3 ./ 3 - 1)) ./ 4 - (nds_2_2 .* (root3 ./ 3 + 1)) ./ 4 + (nds_3_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_4_2 .* (root3 ./ 3 - 1)) ./ 4;
			J(2, 1, :) = (nds_2_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_1_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_3_1 .* (root3 ./ 3 - 1)) ./ 4 + (nds_4_1 .* (root3 ./ 3 - 1)) ./ 4;
			J(2, 2, :) = (nds_2_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_1_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_3_2 .* (root3 ./ 3 - 1)) ./ 4 + (nds_4_2 .* (root3 ./ 3 - 1)) ./ 4;

			%Terms of B3 matrix
			B3(1, 1, :) = root3 ./ 12 - 1 ./ 4;
			B3(1, 4, :) = - root3 ./ 12 - 1 ./ 4;
			B3(1, 7, :) = root3 ./ 12 + 1 ./ 4;
			B3(1, 10, :) = 1 ./ 4 - root3 ./ 12;
			B3(2, 1, :) = - root3 ./ 12 - 1 ./ 4;
			B3(2, 4, :) = root3 ./ 12 + 1 ./ 4;
			B3(2, 7, :) = 1 ./ 4 - root3 ./ 12;
			B3(2, 10, :) = root3 ./ 12 - 1 ./ 4;
			B3(3, 2, :) = root3 ./ 12 - 1 ./ 4;
			B3(3, 5, :) = - root3 ./ 12 - 1 ./ 4;
			B3(3, 8, :) = root3 ./ 12 + 1 ./ 4;
			B3(3, 11, :) = 1 ./ 4 - root3 ./ 12;
			B3(4, 2, :) = - root3 ./ 12 - 1 ./ 4;
			B3(4, 5, :) = root3 ./ 12 + 1 ./ 4;
			B3(4, 8, :) = 1 ./ 4 - root3 ./ 12;
			B3(4, 11, :) = root3 ./ 12 - 1 ./ 4;
			B3(5, 3, :) = root3 ./ 12 - 1 ./ 4;
			B3(5, 6, :) = - root3 ./ 12 - 1 ./ 4;
			B3(5, 9, :) = root3 ./ 12 + 1 ./ 4;
			B3(5, 12, :) = 1 ./ 4 - root3 ./ 12;
			B3(6, 3, :) = - root3 ./ 12 - 1 ./ 4;
			B3(6, 6, :) = root3 ./ 12 + 1 ./ 4;
			B3(6, 9, :) = 1 ./ 4 - root3 ./ 12;
			B3(6, 12, :) = root3 ./ 12 - 1 ./ 4;

			%Terms of N matrix
			N(1, 1, :) = 1 ./ 6;
			N(1, 4, :) = root3 ./ 6 + 1 ./ 3;
			N(1, 7, :) = 1 ./ 6;
			N(1, 10, :) = 1 ./ 3 - root3 ./ 6;
			N(2, 2, :) = 1 ./ 6;
			N(2, 5, :) = root3 ./ 6 + 1 ./ 3;
			N(2, 8, :) = 1 ./ 6;
			N(2, 11, :) = 1 ./ 3 - root3 ./ 6;
			N(3, 3, :) = 1 ./ 6;
			N(3, 6, :) = root3 ./ 6 + 1 ./ 3;
			N(3, 9, :) = 1 ./ 6;
			N(3, 12, :) = 1 ./ 3 - root3 ./ 6;

	%Define matrices that depend on Gauss point
		case 3
			%Terms of Jacobian
			J(1, 1, :) = (nds_1_1 .* (root3 ./ 3 - 1)) ./ 4 - (nds_2_1 .* (root3 ./ 3 + 1)) ./ 4 + (nds_3_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_4_1 .* (root3 ./ 3 - 1)) ./ 4;
			J(1, 2, :) = (nds_1_2 .* (root3 ./ 3 - 1)) ./ 4 - (nds_2_2 .* (root3 ./ 3 + 1)) ./ 4 + (nds_3_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_4_2 .* (root3 ./ 3 - 1)) ./ 4;
			J(2, 1, :) = (nds_1_1 .* (root3 ./ 3 - 1)) ./ 4 - (nds_2_1 .* (root3 ./ 3 - 1)) ./ 4 + (nds_3_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_4_1 .* (root3 ./ 3 + 1)) ./ 4;
			J(2, 2, :) = (nds_1_2 .* (root3 ./ 3 - 1)) ./ 4 - (nds_2_2 .* (root3 ./ 3 - 1)) ./ 4 + (nds_3_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_4_2 .* (root3 ./ 3 + 1)) ./ 4;

			%Terms of B3 matrix
			B3(1, 1, :) = root3 ./ 12 - 1 ./ 4;
			B3(1, 4, :) = - root3 ./ 12 - 1 ./ 4;
			B3(1, 7, :) = root3 ./ 12 + 1 ./ 4;
			B3(1, 10, :) = 1 ./ 4 - root3 ./ 12;
			B3(2, 1, :) = root3 ./ 12 - 1 ./ 4;
			B3(2, 4, :) = 1 ./ 4 - root3 ./ 12;
			B3(2, 7, :) = root3 ./ 12 + 1 ./ 4;
			B3(2, 10, :) = - root3 ./ 12 - 1 ./ 4;
			B3(3, 2, :) = root3 ./ 12 - 1 ./ 4;
			B3(3, 5, :) = - root3 ./ 12 - 1 ./ 4;
			B3(3, 8, :) = root3 ./ 12 + 1 ./ 4;
			B3(3, 11, :) = 1 ./ 4 - root3 ./ 12;
			B3(4, 2, :) = root3 ./ 12 - 1 ./ 4;
			B3(4, 5, :) = 1 ./ 4 - root3 ./ 12;
			B3(4, 8, :) = root3 ./ 12 + 1 ./ 4;
			B3(4, 11, :) = - root3 ./ 12 - 1 ./ 4;
			B3(5, 3, :) = root3 ./ 12 - 1 ./ 4;
			B3(5, 6, :) = - root3 ./ 12 - 1 ./ 4;
			B3(5, 9, :) = root3 ./ 12 + 1 ./ 4;
			B3(5, 12, :) = 1 ./ 4 - root3 ./ 12;
			B3(6, 3, :) = root3 ./ 12 - 1 ./ 4;
			B3(6, 6, :) = 1 ./ 4 - root3 ./ 12;
			B3(6, 9, :) = root3 ./ 12 + 1 ./ 4;
			B3(6, 12, :) = - root3 ./ 12 - 1 ./ 4;

			%Terms of N matrix
			N(1, 1, :) = 1 ./ 3 - root3 ./ 6;
			N(1, 4, :) = 1 ./ 6;
			N(1, 7, :) = root3 ./ 6 + 1 ./ 3;
			N(1, 10, :) = 1 ./ 6;
			N(2, 2, :) = 1 ./ 3 - root3 ./ 6;
			N(2, 5, :) = 1 ./ 6;
			N(2, 8, :) = root3 ./ 6 + 1 ./ 3;
			N(2, 11, :) = 1 ./ 6;
			N(3, 3, :) = 1 ./ 3 - root3 ./ 6;
			N(3, 6, :) = 1 ./ 6;
			N(3, 9, :) = root3 ./ 6 + 1 ./ 3;
			N(3, 12, :) = 1 ./ 6;

	%Define matrices that depend on Gauss point
		case 4
			%Terms of Jacobian
			J(1, 1, :) = (nds_2_1 .* (root3 ./ 3 - 1)) ./ 4 - (nds_1_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_3_1 .* (root3 ./ 3 - 1)) ./ 4 + (nds_4_1 .* (root3 ./ 3 + 1)) ./ 4;
			J(1, 2, :) = (nds_2_2 .* (root3 ./ 3 - 1)) ./ 4 - (nds_1_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_3_2 .* (root3 ./ 3 - 1)) ./ 4 + (nds_4_2 .* (root3 ./ 3 + 1)) ./ 4;
			J(2, 1, :) = (nds_1_1 .* (root3 ./ 3 - 1)) ./ 4 - (nds_2_1 .* (root3 ./ 3 - 1)) ./ 4 + (nds_3_1 .* (root3 ./ 3 + 1)) ./ 4 - (nds_4_1 .* (root3 ./ 3 + 1)) ./ 4;
			J(2, 2, :) = (nds_1_2 .* (root3 ./ 3 - 1)) ./ 4 - (nds_2_2 .* (root3 ./ 3 - 1)) ./ 4 + (nds_3_2 .* (root3 ./ 3 + 1)) ./ 4 - (nds_4_2 .* (root3 ./ 3 + 1)) ./ 4;

			%Terms of B3 matrix
			B3(1, 1, :) = - root3 ./ 12 - 1 ./ 4;
			B3(1, 4, :) = root3 ./ 12 - 1 ./ 4;
			B3(1, 7, :) = 1 ./ 4 - root3 ./ 12;
			B3(1, 10, :) = root3 ./ 12 + 1 ./ 4;
			B3(2, 1, :) = root3 ./ 12 - 1 ./ 4;
			B3(2, 4, :) = 1 ./ 4 - root3 ./ 12;
			B3(2, 7, :) = root3 ./ 12 + 1 ./ 4;
			B3(2, 10, :) = - root3 ./ 12 - 1 ./ 4;
			B3(3, 2, :) = - root3 ./ 12 - 1 ./ 4;
			B3(3, 5, :) = root3 ./ 12 - 1 ./ 4;
			B3(3, 8, :) = 1 ./ 4 - root3 ./ 12;
			B3(3, 11, :) = root3 ./ 12 + 1 ./ 4;
			B3(4, 2, :) = root3 ./ 12 - 1 ./ 4;
			B3(4, 5, :) = 1 ./ 4 - root3 ./ 12;
			B3(4, 8, :) = root3 ./ 12 + 1 ./ 4;
			B3(4, 11, :) = - root3 ./ 12 - 1 ./ 4;
			B3(5, 3, :) = - root3 ./ 12 - 1 ./ 4;
			B3(5, 6, :) = root3 ./ 12 - 1 ./ 4;
			B3(5, 9, :) = 1 ./ 4 - root3 ./ 12;
			B3(5, 12, :) = root3 ./ 12 + 1 ./ 4;
			B3(6, 3, :) = root3 ./ 12 - 1 ./ 4;
			B3(6, 6, :) = 1 ./ 4 - root3 ./ 12;
			B3(6, 9, :) = root3 ./ 12 + 1 ./ 4;
			B3(6, 12, :) = - root3 ./ 12 - 1 ./ 4;

			%Terms of N matrix
			N(1, 1, :) = 1 ./ 6;
			N(1, 4, :) = 1 ./ 3 - root3 ./ 6;
			N(1, 7, :) = 1 ./ 6;
			N(1, 10, :) = root3 ./ 6 + 1 ./ 3;
			N(2, 2, :) = 1 ./ 6;
			N(2, 5, :) = 1 ./ 3 - root3 ./ 6;
			N(2, 8, :) = 1 ./ 6;
			N(2, 11, :) = root3 ./ 6 + 1 ./ 3;
			N(3, 3, :) = 1 ./ 6;
			N(3, 6, :) = 1 ./ 3 - root3 ./ 6;
			N(3, 9, :) = 1 ./ 6;
			N(3, 12, :) = root3 ./ 6 + 1 ./ 3;

	end

	%Determinant of Jacobian
	detJ = J(1, 1, :) .* J(2, 2, :) - J(1, 2, :) .* J(2, 1, :);

	%Terms of B2 matrix
	B2(1, 1, :) = J(2, 2, :) ./ detJ;
	B2(1, 2, :) = -J(1, 2, :) ./ detJ;
	B2(2, 1, :) = -J(2, 1, :) ./ detJ;
	B2(2, 2, :) = J(1, 1, :) ./ detJ;
	B2(4, 3, :) = J(2, 2, :) ./ detJ;
	B2(4, 4, :) = -J(1, 2, :) ./ detJ;
	B2(5, 3, :) = -J(2, 1, :) ./ detJ;
	B2(5, 4, :) = J(1, 1, :) ./ detJ;
	B2(7, 5, :) = J(2, 2, :) ./ detJ;
	B2(7, 6, :) = -J(1, 2, :) ./ detJ;
	B2(8, 5, :) = -J(2, 1, :) ./ detJ;
	B2(8, 6, :) = J(1, 1, :) ./ detJ;

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
