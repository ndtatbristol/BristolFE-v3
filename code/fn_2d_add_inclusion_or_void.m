function mod = fn_2d_add_inclusion_or_void(mod, matls, el_types, scat_pts, scat_matl)
%USAGE
%   mod = fn_2d_add_inclusion_or_void(mod, matls, el_types, scat_pts, scat_matl)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Adds scatterer to existing model by turning all elements inside
%   scat_pts to either matl(scat_matl) or void if = scat_matl
%INPUTS
%   mod - existing model
%   matls - cell array of materials in model
%   el_types - cell array of element types in model
%   scat_pts - n_ptsx2 matrix of points describing a closed boundary around
%   scatterer or void
%   scat_matl - material index of material to go inside scatterer region or
%   zero to create a void
%OUTPUT
%   mod - modified model with scatterer
%NOTES
%   This looks like it was written to handle sub-domains (extra variables
%   involving node numbers are updated) and fluid-solid interfaces. Not
%   clear if either of these are actually required and it would be neater
%   if they were not. This should be checked.
%   7/10/25 - removed lines that remove and then add fluid-solid interface
%   elements on either side and it still seems ok. For subdoms, it would be
%   neater to just return old_nds as second parameter than can be used if
%   needed to update the extra fields needed for subdomain models.
%--------------------------------------------------------------------------


interface_el_name = 'ASI2D2';

%Remove interface elements if there are any
%mod = fn_remove_fluid_solid_interface_els(mod, el_types);
% els_in_use = ~strcmp(mod.el_typ_i, interface_el_name);
% [~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);


if scat_matl > 0
    mod = rmfield(mod, 'el_typ_i');
    mod = fn_set_els_inside_bdry_to_mat(mod, scat_pts, scat_matl);
else
    % [~, els_in_use] = fn_elements_in_region2(mod.nds, mod.els, scat_pts);
    [~, els_in_use] = fn_elements_in_region(mod, scat_pts);
    [~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);
    [mod.nds, mod.els, old_nds] = fn_remove_unused_nodes(mod.nds, mod.els);
    %Following needed for sub-domain models?
    if isfield(mod, 'bdry_lyrs')
        mod.bdry_lyrs = mod.bdry_lyrs(old_nds);
    end
    if isfield(mod, 'main_nd_i')
        mod.main_nd_i = mod.main_nd_i(old_nds);
    end
end

%Add interface elements if needed
%mod = fn_add_fluid_solid_interface_els(mod, el_types);

if isfield(mod, 'inner_bndry_pts')
    %Set flag on which elements are within domain
    mod.int_el_i = fn_elements_in_region(mod, mod.inner_bndry_pts);
end

end