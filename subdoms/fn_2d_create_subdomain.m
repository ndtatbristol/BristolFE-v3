function dm_mod = fn_2d_create_subdomain(mn_mod, inner_bndry_pts, abs_layer_thick)
%USAGE
%   dm_mod = fn_2d_create_subdomain(mn_mod, inner_bdry_vtcs, abs_layer_thick)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Creates a 2D subdomain model from a larger 2D model with the subdomain
%   defined by a series of points describing its boundary in the main model
%   and the thickness of absorbing layer to use around that boundary in the 
%   subdomain model.
%INPUTS
%   mn_mod - structure describing the main model which must include the
%   following fields:
%       TODO
%   inner_bndry_pts - [n_inner_bndry_pts  x 2] matrix of points
%   describing inner boundary of subdomain
%   abs_layer_thick - thickness of the absorbing layer around the subdomain
%OUTPUTS
%   dm_mod - structure describing the subdomain model with fields
%       nds - [n_nds x 2] matrix of coordinates of subdomain model nodes
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
%       inner_bndry_pts - [n_inner_bndry_pts  x 2] pass-through of inner_bndry_pts
% %       outer_bndry_vtcs - [n_outer_bdry_pts  x 2] vertices describing
% %       outer boundary of subdomain model
%--------------------------------------------------------------------------


dm_mod = fn_create_subdomain(mn_mod, inner_bndry_pts, [], abs_layer_thick);
dm_mod.inner_bndry_pts = inner_bndry_pts;

end
