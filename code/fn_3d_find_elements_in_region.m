function [in, out] = fn_3d_find_elements_in_region(mod, vtcs, fcs, varargin)
%USAGE
%   [in, out] = fn_3d_find_elements_in_region(mod, vtcs, fcs [, interior_pt])
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
%   [interior_pt - coordinates of a point that is inside region. If not
%   specified, it is assumed that interior_pt = mean(vtcs) which may not be
%   correct for oddly-shaped regions, e.g. toroids.

%TODO - interior point should not be necessary on any of these functions -
%just pick a test point that is outside region bounded by vertices and make
%sure that is outside. Point to apply this is in the
%fn_signed_distance_to_boundary function to force interior points to be
%negative distances.

%OUTPUTS
%   in - [n_els x 1] matrix of binary values indicating whether each
%   element (based on its centre) is inside of outside the region
%   out - [n_els x 1] = ~in.
%--------------------------------------------------------------------------
if numel(varargin) >= 1
    interior_pt = varargin{1};
else
    interior_pt = mean(vtcs);
end

el_centres = fn_calc_element_centres(mod.nds, mod.els);
%reduce search space by excluding all elements outside bounding box of
%vertices
min_vtcs = min(vtcs);
max_vtcs = max(vtcs);
possible_els_in_region = all((el_centres <= max_vtcs) & (el_centres >= min_vtcs), 2);

d = zeros(size(el_centres, 1), 1);
d(possible_els_in_region) = fn_signed_dist_to_bdry(el_centres(possible_els_in_region, :), vtcs, fcs, interior_pt);
in = d < 0;
out = ~in;
end