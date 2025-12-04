function dm_mod = fn_3d_create_subdomain(mn_mod, el_types, inner_bndry_vtcs, inner_bndry_fcs, abs_layer_thick)
%USAGE
%   dm_mod = fn_3d_create_subdomain(mn_mod, inner_bndry_vtcs, inner_bndry_fcs, abs_layer_thick)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Creates a 3D subdomain model from a larger 3D model with the subdomain
%   defined by a series of vertices and faces describing its boundary in 
%   the main model and the thickness of absorbing layer to use around that 
%   boundary in the subdomain model.
%INPUTS
%   mn_mod - structure describing the main model which must include the
%   following fields:
%       TODO
%   inner_bndry_vtcs - [n_inner_bndry_vtcs  x 3] matrix of vertices
%   describing inner boundary of subdomain
%   inner_bndry_fcs - [n_inner_bndry_fcs  x 3] matrix of vertex indices
%   describing inner boundary faces
%   abs_layer_thick - thickness of the absorbing layer around the subdomain
%OUTPUTS
%   dm_mod - structure describing the subdomain model with fields
%       nds - [n_nds x 3] matrix of coordinates of subdomain model nodes
%       els - [n_els x max_nds_per_el] matrix of nodes for each element in
%       subdomain model
%       el_mat_i, el_typ_i, el_abs_i - [n_els x 1] vectors of element
%       material indices, element type indices, and element absorbing
%       indices for each element in subdomain model
%       bdry_lyrs - [n_nds x 1] vector containing 1,2,3, or 4 for nodes in
%       boundary layers (1 is innermost) of subdomain model or 0 for nodes 
%       that are not
%       main_nd_i - [n_nds x 1] vector of equivalent node numbers in main
%       model for each node in subdomain model
%       inner_bndry_vtcs - [n_inner_bndry_vtcs x 3] pass-through of inner_bndry_pts
%       inner_bndry_fcs - [n_inner_bndry_fcs x 3] pass-through of inner_bndry_fcs
% %       outer_bndry_vtcs - [n_outer_bdry_pts  x 2] vertices describing
% %       outer boundary of subdomain model
%--------------------------------------------------------------------------


dm_mod = fn_create_subdomain(mn_mod, el_types, inner_bndry_vtcs, inner_bndry_fcs, abs_layer_thick);
dm_mod.inner_bndry_pts = inner_bndry_vtcs;
dm_mod.inner_bndry_fcs = inner_bndry_fcs;
end

