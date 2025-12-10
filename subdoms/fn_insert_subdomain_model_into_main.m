function [val_mod, old_nds, new_nds] = fn_insert_subdomain_model_into_main(mn_mod, dm_mod)

val_mod = mn_mod;
ndim = size(val_mod.nds, 2);

%Remove the elements from main model that are inside region
switch ndim
    case 2
        els_in_use = ~fn_2d_find_elements_in_region(val_mod, dm_mod.inner_bndry_pts);
    case 3
        els_in_use = ~fn_3d_find_elements_in_region(val_mod, dm_mod.inner_bndry_pts, dm_mod.inner_bndry_fcs);
end
[~, ~, val_mod.els, val_mod.el_mat_i, val_mod.el_abs_i, val_mod.el_typ_i] = fn_remove_unused_elements(els_in_use, val_mod.els, val_mod.el_mat_i, val_mod.el_abs_i, val_mod.el_typ_i);

%Note that above step removes elements from main model, but not nodes
%(those within the subdomain region will not be associated with any
%elements. This is important to preserve the node numbering in the main
%model as the subdomain model references this through the main_nd_i
%parameter.

%remove elements in subdomain that are outside region (i.e. the absorbing
%boundaries)
switch ndim
    case 2
        els_in_use = fn_2d_find_elements_in_region(dm_mod, dm_mod.inner_bndry_pts);
    case 3
        els_in_use = fn_3d_find_elements_in_region(dm_mod, dm_mod.inner_bndry_pts, dm_mod.inner_bndry_fcs);
end
[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i] = fn_remove_unused_elements(els_in_use, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i);

%Now the sub-domain also has a load of unattached nodes that were in the
%absorbing boundaries. 

%Stick the subdomain nodes at the end of the main model nodes and remember
%the offset that needs to be added to references to these
dm_nd_offset = size(val_mod.nds, 1);
val_mod.nds = [val_mod.nds; dm_mod.nds];

%Find the boundary nodes (dm_mod.bdry_lyrs == 1, which correspond to remaining nodes that
%are common to main and subdomain models). It would feel better to do this
%by checking for physical overlap since from above logic we have guarantee
%of contiguous elements in modified main and subdom models, but no
%guarantee from within this function that boundary corresponds to dm_mod.bdry_lyrs == 1
%... that comes from logic of how subdom is defined in the first place.
bdry_nds = find(dm_mod.bdry_lyrs == 1);
i = ismember(dm_mod.els, bdry_nds);

%change interface nodes in subdomain model elements to main model nodes and add
%offset onto the other nodes in subdomain model elements. Keep track of
%zeros as these need to remain zeros afterwards (and not have the offset
%added on)
j = dm_mod.els == 0;
dm_mod.els(i) = fn_remap_matrix(dm_mod.els(i), dm_mod.main_nd_i);
dm_mod.els(~i) = dm_mod.els(~i) + dm_nd_offset;
dm_mod.els(j) = 0;

%append the domain model els to the main model, adding dm_nd_offset to the node indices
%as the subdomain model nodes were appended to the main model nodes
val_mod.els =      [val_mod.els;      dm_mod.els];
val_mod.el_mat_i = [val_mod.el_mat_i; dm_mod.el_mat_i];
val_mod.el_abs_i = [val_mod.el_abs_i; dm_mod.el_abs_i];
val_mod.el_typ_i = [val_mod.el_typ_i; dm_mod.el_typ_i];

%Drop the unused nodes and update elements
[val_mod.nds, val_mod.els, ~, new_nds1] = fn_remove_unused_nodes(val_mod.nds, val_mod.els);

%Sort nodes and update elements
[val_mod.nds, tmp] = sortrows(val_mod.nds, (size(val_mod.nds,2):-1:1));
new_nds2 = fn_inverse_map(tmp);
val_mod.els = fn_remap_matrix(val_mod.els, new_nds2);

%Work out the old_nds and new_nds return values
j = new_nds1 == 0;
new_nds1(j) = 1;
new_nds = new_nds2(new_nds1);
old_nds = fn_inverse_map(new_nds);
new_nds(j) = 0;
end


