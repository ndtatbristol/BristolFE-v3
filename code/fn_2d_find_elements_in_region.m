function [in, out] = fn_2d_find_elements_in_region(mod, region)
%USAGE
%   [in, out] = fn_2d_find_elements_in_region(mod, region)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Returns logical n_els x 1 vectors indicating whether elements in model
%   are inside or outside the specified region.
%INPUTS
%   mod - structured variable describing 3D model geometry which must
%   contain fields:
%       nds - [n_nds x 2] matrix of nodal coordinates
%       els - [n_els x max_nds_per_el] matrix of node indices for each
%       element in model
%   region - [n_vtcs x 2] matrix of coordinates of polygon defining region
%OUTPUTS
%   in - [n_els x 1] matrix of binary values indicating whether each
%   element (based on its centre) is inside of outside the region
%   out - [n_els x 1] = ~in.
%--------------------------------------------------------------------------

el_centres = fn_calc_element_centres(mod.nds, mod.els);
in = inpolygon(el_centres(:,1), el_centres(:,2), region(:,1), region(:,2));
out = ~in;
end