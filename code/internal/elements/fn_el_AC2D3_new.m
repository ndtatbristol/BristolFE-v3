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
%	Paul Wilcox (27-Mar-2026 20:42:59)

root3 = sqrt(3);
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
	el_Q = [];
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


%Define Gauss points and weights
Q = [0.3333333333 0.3333333333 ];
W = [0.5];

%Zero the outputs
el_K = zeros(size(els, 1), 3, 3);
tmp_M = zeros(size(els, 1), 3, 3);
el_C = zeros(size(els, 1), 3, 3);
BTDB = zeros(size(els, 1), 3, 3);
BTD = zeros(size(els, 1), 3, 1);
NTN = zeros(size(els, 1), 3, 3);

%Loop over Gauss points
for i = 1:size(Q, 1)
    q1 = Q(i, 1);
    q2 = Q(i, 2);

    %Evaluate Jacobian
    detJ = nds_1_1 .* nds_2_2 - nds_1_2 .* nds_2_1 - nds_1_1 .* nds_3_2 + nds_1_2 .* nds_3_1 + nds_2_1 .* nds_3_2 - nds_2_2 .* nds_3_1;

    %Evaluate B matrix
    B = zeros(size(els, 1), 1, 3);
    B(:, 1, 1) = (-((nds_1_1 - nds_2_1) ./ detJ - (nds_1_2 - nds_2_2) ./ detJ - (nds_1_1 - nds_3_1) ./ detJ + (nds_1_2 - nds_3_2) ./ detJ) ./ (D .* rho)) .^ (1 ./ 2);
    B(:, 1, 2) = (-((nds_1_1 - nds_3_1) ./ detJ - (nds_1_2 - nds_3_2) ./ detJ) ./ (D .* rho)) .^ (1 ./ 2);
    B(:, 1, 3) = (((nds_1_1 - nds_2_1) ./ detJ - (nds_1_2 - nds_2_2) ./ detJ) ./ (D .* rho)) .^ (1 ./ 2);


    %Evaluate B'D
    BTD(:, 1, 1) =  + B(:, 1, 1) * D(1, 1);
    BTD(:, 2, 1) =  + B(:, 1, 2) * D(1, 1);
    BTD(:, 3, 1) =  + B(:, 1, 3) * D(1, 1);

    %Evaluate B'DB
    BTDB(:, 1, 1) =  + BTD(:, 1, 1) .* B(:, 1, 1);
    BTDB(:, 1, 2) =  + BTD(:, 1, 1) .* B(:, 1, 2);
    BTDB(:, 1, 3) =  + BTD(:, 1, 1) .* B(:, 1, 3);
    BTDB(:, 2, 2) =  + BTD(:, 2, 1) .* B(:, 1, 2);
    BTDB(:, 2, 3) =  + BTD(:, 2, 1) .* B(:, 1, 3);
    BTDB(:, 3, 3) =  + BTD(:, 3, 1) .* B(:, 1, 3);

    %Evaluate contribution to K at Gauss point and accumulate
    el_K(:, 1, 1) = el_K(:, 1, 1) + BTDB(:, 1, 1) .* detJ * W(i);
    el_K(:, 1, 2) = el_K(:, 1, 2) + BTDB(:, 1, 2) .* detJ * W(i);
    el_K(:, 1, 3) = el_K(:, 1, 3) + BTDB(:, 1, 3) .* detJ * W(i);
    el_K(:, 2, 2) = el_K(:, 2, 2) + BTDB(:, 2, 2) .* detJ * W(i);
    el_K(:, 2, 3) = el_K(:, 2, 3) + BTDB(:, 2, 3) .* detJ * W(i);
    el_K(:, 3, 3) = el_K(:, 3, 3) + BTDB(:, 3, 3) .* detJ * W(i);

    %Evaluate N matrix
    N = zeros(size(els, 1), 1, 3);
    N(:, 1, 1) = ((q1 + q2 - 1) ./ (D .* rho)) .^ (1 ./ 2);
    N(:, 1, 2) = (-q1 ./ (D .* rho)) .^ (1 ./ 2);
    N(:, 1, 3) = (-q2 ./ (D .* rho)) .^ (1 ./ 2);


    %Evaluate N'N
    NTN(:, 1, 1) =  + N(:, 1, 1) .* N(:, 1, 1);
    NTN(:, 1, 2) =  + N(:, 1, 1) .* N(:, 1, 2);
    NTN(:, 1, 3) =  + N(:, 1, 1) .* N(:, 1, 3);
    NTN(:, 2, 2) =  + N(:, 1, 2) .* N(:, 1, 2);
    NTN(:, 2, 3) =  + N(:, 1, 2) .* N(:, 1, 3);
    NTN(:, 3, 3) =  + N(:, 1, 3) .* N(:, 1, 3);

    %Evaluate contribution to M at Gauss point and accumulate
    tmp_M(:, 1, 1) = tmp_M(:, 1, 1) + NTN(:, 1, 1) .* rho * W(i);
    tmp_M(:, 1, 2) = tmp_M(:, 1, 2) + NTN(:, 1, 2) .* rho * W(i);
    tmp_M(:, 1, 3) = tmp_M(:, 1, 3) + NTN(:, 1, 3) .* rho * W(i);
    tmp_M(:, 2, 2) = tmp_M(:, 2, 2) + NTN(:, 2, 2) .* rho * W(i);
    tmp_M(:, 2, 3) = tmp_M(:, 2, 3) + NTN(:, 2, 3) .* rho * W(i);
    tmp_M(:, 3, 3) = tmp_M(:, 3, 3) + NTN(:, 3, 3) .* rho * W(i);
end


%Copy upper triangles of K and tmp_M into lower for symmmetry
el_K(:, 2, 1) = el_K(:, 1, 2);
tmp_M(:, 2, 1) = tmp_M(:, 1, 2);
el_K(:, 3, 1) = el_K(:, 1, 3);
tmp_M(:, 3, 1) = tmp_M(:, 1, 3);
el_K(:, 3, 2) = el_K(:, 2, 3);
tmp_M(:, 3, 2) = tmp_M(:, 2, 3);

%Diagonalise M
M = zeros(size(els, 1), 3, 3);
el_M(:, 1, 1) = sum(tmp_M(:, :, 1));
el_M(:, 2, 2) = sum(tmp_M(:, :, 2));
el_M(:, 3, 3) = sum(tmp_M(:, :, 3));
%CRemove unwanted DOFs from element matrices
[loc_nd, loc_df, el_K, el_C, el_M] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use, el_K, el_C, el_M);

end
