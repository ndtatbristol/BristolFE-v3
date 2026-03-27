function [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_CPE4_test3(nds, els, D, rho, varargin)
%SUMMARY
%	This function was created automatically by fn_create_element_matrix_file
%	and contains code to return the stiffness and mass matrices
%	for multiple elements of the same material and type given by the latter
%	part of the filename, fn_el_CPE4_test3.
%INPUTS
%	nds - n_nds x n_dims matrix of nodal coordinates
%	els - n_els x n_nds_per_el matrix of node indices for each elements
%	D - ns x ns material stiffness matrix
%	rho - material density
%	[dofs_to_use = [] - optional string listing the DoFs to use, e.g. '12'. Use [] for all]
%OUTPUTS
%	el_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices
%AUTHOR
%	Paul Wilcox (26-Mar-2026 08:45:27)

%Deal with optional argument about which DOFs to use
if isempty(varargin)
	dofs_to_use = [];
else
	dofs_to_use = varargin{1};
end

%Record the local node numbers of the element stiffness matrices
