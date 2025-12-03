function mod = fn_3d_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i, varargin)
%USAGE
%   mod = fn_3d_add_inclusion_or_void(mod, el_types, scat_pts, scat_matl, scat_el_typ[, scat_interior_pt])
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Adds scatterer to existing model by turning all elements inside
%   scat_pts to either matl(scat_matl) or void if = scat_matl
%INPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, and element nodes, mod.els.
%   el_types - cell array of element type names to be used in model
%   scat_vtcs - n_ptsx3 matrix of points describing vertices of scatterer
%   scat_fcs - m_fcsx3 matrix of vertex indices of triangular facets
%   describing scatterer
%   scat_matl_i - material index of material to go inside scatterer region or
%   zero to create a void
%   scat_el_typ_i - element type index of material inside scatterer region
%   (ignored if scat_matl == 0)
%   [scat_interior_pt - a point inside the scatterer to define the interior
%   of 3D scatterers. If not specified, the mean scatterer coordinate will 
%   be assumed to be inside ... but this may not be the case for all
%   scatterer shapes, e.g. donuts!]
%OUTPUT
%   mod - modified model with scatterer
%NOTES
%   This handles subdomains by checking for existance of other fields in
%   mod that involve node numbers as these have to be updated if nodes are
%   removed (i.e. when a void is added; an inclusion does not remove any
%   nodes or elements, it just changes element material)
%--------------------------------------------------------------------------
if numel(varargin) >= 1
    scat_interior_pt = varargin{1};
else
    scat_interior_pt = mean(mod.nds);
end

mod = fn_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i, scat_interior_pt);


end