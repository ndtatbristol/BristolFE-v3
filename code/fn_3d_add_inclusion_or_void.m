function mod = fn_3d_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i)
%USAGE
%   mod = fn_3d_add_inclusion_or_void(mod, el_types, scat_pts, scat_matl, scat_el_typ)
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
%OUTPUT
%   mod - modified model with scatterer
%NOTES
%   This handles subdomains by checking for existance of other fields in
%   mod that involve node numbers as these have to be updated if nodes are
%   removed (i.e. when a void is added; an inclusion does not remove any
%   nodes or elements, it just changes element material)
%--------------------------------------------------------------------------

mod = fn_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i);


end