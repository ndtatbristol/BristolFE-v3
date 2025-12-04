function mod = fn_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i)
%USAGE
%   mod = fn_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl, scat_el_typ)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Generic function for 2D or 3D
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
ndims = size(mod.nds, 2);
if ndims == 2
    els_in_inclusion = fn_2d_find_elements_in_region(mod, scat_vtcs);
else
    els_in_inclusion = fn_3d_find_elements_in_region(mod, scat_vtcs, scat_fcs);
end
if scat_matl_i > 0
    %only change material and type of non-interface elements, otherwise
    %interface elements cannot be detected and deleted later
    int_el_typ_i = fn_el_typ_indices_for_class(el_types, 'fluid_solid_interface');
    for i = 1:numel(int_el_typ_i)
        els_in_inclusion = els_in_inclusion & (mod.el_typ_i ~= int_el_typ_i(i));
    end
    mod.el_mat_i(els_in_inclusion) = scat_matl_i;
    mod.el_typ_i(els_in_inclusion) = scat_el_typ_i;
    %Add interface elements (may not be necessary always, but no harm in
    %calling as they will be necessary if it is a solid inclusion in a
    %liquid or vice versa
    mod = fn_add_fluid_solid_interface_els(mod, el_types);
else
    [~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(~els_in_inclusion, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);
    [mod.nds, mod.els, old_nds] = fn_remove_unused_nodes(mod.nds, mod.els);
    %Following needed for sub-domain models
    if isfield(mod, 'bdry_lyrs')
        mod.bdry_lyrs = mod.bdry_lyrs(old_nds);
    end
    if isfield(mod, 'main_nd_i')
        mod.main_nd_i = mod.main_nd_i(old_nds);
    end
end

% %Following needed for sub-domain models
% if isfield(mod, 'inner_bndry_pts')
%     %Set flag on which elements are within domain
%     mod.int_el_i = fn_2d_find_elements_in_region(mod, mod.inner_bndry_pts);
% end

end