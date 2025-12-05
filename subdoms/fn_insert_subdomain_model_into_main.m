function mod = fn_insert_subdomain_model_into_main(mn_mod, dm_mod)

mod = mn_mod;
ndim = size(mod.nds, 2);

%Remove the elements from main model that are inside region
switch ndim
    case 2
        els_in_use = ~fn_2d_find_elements_in_region(mod, dm_mod.inner_bndry_pts);
    case 3
        els_in_use = ~fn_3d_find_elements_in_region(mod, dm_mod.inner_bndry_pts, dm_mod.inner_bndry_fcs);
end
[~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);

%Stick the subdomain nodes at the end of the main model nodes and remember
%the offset that needs to be added to references to these
dm_nd_offset = size(mod.nds, 1);
mod.nds = [mod.nds; dm_mod.nds];

%remove elements in subdomain that are outside region
% els_in_use = fn_elements_in_region(dm_mod, dm_mod.inner_bndry_pts);
% els_in_use = fn_2d_find_elements_in_region(dm_mod, dm_mod.inner_bndry_pts);
switch ndim
    case 2
        els_in_use = fn_2d_find_elements_in_region(dm_mod, dm_mod.inner_bndry_pts);
    case 3
        els_in_use = fn_3d_find_elements_in_region(dm_mod, dm_mod.inner_bndry_pts, dm_mod.inner_bndry_fcs);
end
[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i] = fn_remove_unused_elements(els_in_use, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i);

%Change boundary node references in sub-domain model to main node
%equivalents
ed_nds = find(dm_mod.bdry_lyrs == 1);
i = ismember(dm_mod.els, ed_nds);
dm_mod.els(i) = interp1(1:size(dm_mod.nds, 1), dm_mod.main_nd_i, dm_mod.els(i), 'nearest') - dm_nd_offset;
dm_mod.els(dm_mod.els == 0) = -dm_nd_offset; %need to handle zeros (which need to be still zeros when offset)
mod.els =      [mod.els;      dm_mod.els + dm_nd_offset];
mod.el_mat_i = [mod.el_mat_i; dm_mod.el_mat_i];
mod.el_abs_i = [mod.el_abs_i; dm_mod.el_abs_i];
mod.el_typ_i = [mod.el_typ_i; dm_mod.el_typ_i];

end


