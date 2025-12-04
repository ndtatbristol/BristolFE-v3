function [in, out] = fn_3d_find_elements_in_region(mod, vtcs, fcs)
%USAGE
%   [in, out] = fn_3d_find_elements_in_region(mod, vtcs, fcs)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Returns logical n_els x 1 vectors indicating whether elements in model
%   are inside or outside the specified region, which is defined in terms 
%   of vertices and triangular faces.
%INPUTS
%   mod - structured variable describing 3D model geometry which must
%   contain fields:
%       nds - [n_nds x 3] matrix of nodal coordinates
%       els - [n_els x max_nds_per_el] matrix of node indices for each
%       element in model
%   vtcs - [n_vtcs x 3] matrix of coordinates of vertices defining region
%   fcs - [n_fcs x 3] matrix of vertex indices defining faces of region
%OUTPUTS
%   in - [n_els x 1] matrix of binary values indicating whether each
%   element (based on its centre) is inside of outside the region
%   out - [n_els x 1] = ~in.
%--------------------------------------------------------------------------
el_centres = fn_calc_element_centres(mod.nds, mod.els);

%reduce search space by excluding all elements outside bounding box of
%vertices
min_vtcs = min(vtcs);
max_vtcs = max(vtcs);
possible_els_in_region = all((el_centres <= max_vtcs) & (el_centres >= min_vtcs), 2);

d = zeros(size(el_centres, 1), 1);
d(possible_els_in_region) = fn_3d_signed_dist_to_bdry(el_centres(possible_els_in_region, :), vtcs, fcs);
in = d < 0;
out = ~in;
end