function [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_AC2D3_new(nds, els, D, rho, varargin)
%SUMMARY
%	This function was created automatically by fn_create_element_matrix_file
%	and contains code to return the stiffness and mass matrices
%	for multiple elements of the same material and type given by the latter
%	part of the filename, fn_el_AC2D3_new.
%INPUTS
%	nds - n_nds x n_dims matrix of nodal coordinates
%	els - n_els x n_nds_per_el matrix of node indices for each elements
%	D - ns x ns material stiffness matrix
%	rho - material density
%	[dofs_to_use = [] - optional string listing the DoFs to use, e.g. '12'. Use [] for all]
%OUTPUTS
%	el_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices
%AUTHOR
%	Paul Wilcox (29-Mar-2026 21:38:27)

rt3 = sqrt(3);
%Deal with optional argument about which DOFs to use
if isempty(varargin)
	dofs_to_use = [];
else
	dofs_to_use = varargin{1};
end

%Record the local node numbers of the element stiffness matrices
loc_nd = [1  2  3];

%Record the local DOFs of the element stiffness matrices
loc_df = [4  4  4];

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

%Temporary matrices of nodal coordinates to save time
nds_1_1 = nds(els(:, 1), 1);
nds_1_2 = nds(els(:, 1), 2);
nds_2_1 = nds(els(:, 2), 1);
nds_2_2 = nds(els(:, 2), 2);
nds_3_1 = nds(els(:, 3), 1);
nds_3_2 = nds(els(:, 3), 2);


%Some constants
no_gauss_pts = 1;
no_els = size(els, 1);


%Zero the outputs
el_K = zeros(3, 3, no_els);
el_M_tmp = zeros(3, 3, no_els);
el_C = zeros(3, 3, no_els);

detJ = zeros(no_els, 1);
B = zeros(1, 3, no_els);
N = zeros(1, 3, no_els);
%Loop over Gauss points
for i = 1:no_gauss_pts

    %Jacobians, N- and B-matrices at each Gauss point

    switch i
        case 1
            detJ = (nds_1_1 .* nds_2_2) ./ 2 - (nds_1_2 .* nds_2_1) ./ 2 - (nds_1_1 .* nds_3_2) ./ 2 + (nds_1_2 .* nds_3_1) ./ 2 + (nds_2_1 .* nds_3_2) ./ 2 - (nds_2_2 .* nds_3_1) ./ 2;

