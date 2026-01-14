function mod = fn_3d_add_absorbing_layer(mod, abs_bdry_nds, abs_bdry_fcs, abs_bdry_thickness, varargin)
%USAGE
%   mod = fn_3d_add_absorbing_layer(mod, abs_bdry_nds, abs_bdry_fcs, abs_bdry_thickness [, els_to_apply_to])
%AUTHOR
%   Paul Wilcox (2026)
%SUMMARY
%   Adds an absorbing boundary by increasing element absorbing indices
%   proportional to their distance from the specified boundary divided by
%   the specified absorbing boundary layer thickness (i.e. so it reaches
%   one when the distance is equal to the absorbing boundary layer 
%   thickness. The boundary defines the start of the absorbing layer;
%   within the boundary the absorbing index is set to zero.
%INPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, and element nodes, mod.els
%   abs_bdry_pts - coordinates of points defining start absorbing boundary.
%   Does not have to be contained fully within model if absorbing
%   boundaries are only required on certain sides but not others
%   abs_bdry_thickness - absorption index for elements outside boundary
%   will be set to their distance from boundary divided by this value (and
%   capped at unity). Typically this is the distance of the abs_bdry_pts
%   from the edge of the model.
%   [els_to_apply_to - logical vector of size n_els x 1 which specifies
%   which elements in model are affected. Default is all of them, but this
%   parameter enables absorbing boundary to be just applied to part of it,
%   which may be useful in complex geometries]
%OUTPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, element nodes, mod.els, and element absorption
%   index, mod.el_abs_i. 
%NOTES
%   mod.el_abs_i(els_to_apply_to) will be overwritten by this function. 
%   Therefore, if called without els_to_apply_to specified it will be applied
%   to the whole model and will overwrite any previously defined absorption 
%   values.
%--------------------------------------------------------------------------

if numel(varargin) < 1
    els_to_apply_to = ones(size(mod.el_abs_i));
else
    els_to_apply_to = varargin{1};
end

els_to_apply_to = logical(els_to_apply_to);

el_ctrs = fn_calc_element_centres(mod.nds, mod.els(els_to_apply_to, :));
d = fn_3d_signed_dist_to_bdry(el_ctrs, abs_bdry_nds, abs_bdry_fcs);

mod.el_abs_i(els_to_apply_to) = d / abs_bdry_thickness;
mod.el_abs_i(els_to_apply_to & (mod.el_abs_i < 0)) = 0;
mod.el_abs_i(els_to_apply_to & (mod.el_abs_i > 1)) = 1;

end