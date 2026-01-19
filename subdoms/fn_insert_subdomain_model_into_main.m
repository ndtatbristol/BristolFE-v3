function [val_mod, old_nds, new_nds, dm_mod] = fn_insert_subdomain_model_into_main(mn_mod, dm_mod)

%18/1/26 - idea for better way of doing this, that avoids having to use
%fn_{2d/3d}_find_elements_in_region again:
%Use fact that dm_mod.main_nd_i relates nodes dm_mod to mn_mod AND that
%that dm_mod.bdry_lyrs == 1 identifies nodes in dm_mod which define the
%stiching line. Logic something like:
%Define N = find(dm_mod.bdry_lyrs == 1)
%Sub-domain model - identify all els in contact with N and then subset of 
%these in contact with dm_mod.bdry_lyrs == 2.
%Iteratively work outwards to flag (and remove) all elements outside N
%Main model - same idea but work inwards from equivalent nodes to N to
%remove elements inside N.

val_mod = fn_trim_subdomain_from_main(mn_mod, dm_mod);
dm_mod = fn_trim_out_from_subdomain(dm_mod);

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

j = new_nds1 == 0;
new_nds1(j) = 1;
old_nds = fn_inverse_map(new_nds1);
new_nds = new_nds1;
new_nds(j) = 0;
end

function dm_mod = fn_trim_out_from_subdomain(dm_mod)
bdry_nds1 = find(dm_mod.bdry_lyrs == 1);
tmp_els1 = fn_els_on_nds(dm_mod.els, bdry_nds1);
nds_in = find(dm_mod.bdry_lyrs == 2);
tmp_els2 = fn_els_on_nds(dm_mod.els, nds_in);
els_in = tmp_els1 & tmp_els2;
all_els_in_subdomain_abs_bdry = fn_all_els_inside_bdry(dm_mod.els, nds_in, els_in);
[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i] = fn_remove_unused_elements(~all_els_in_subdomain_abs_bdry, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i);
end

function mn_mod = fn_trim_subdomain_from_main(mn_mod, dm_mod)
bdry_nds1 = dm_mod.main_nd_i(dm_mod.bdry_lyrs == 1);
tmp_els1 = fn_els_on_nds(mn_mod.els, bdry_nds1);
bdry_nds2 = dm_mod.main_nd_i(dm_mod.bdry_lyrs == 2);
tmp_els2 = fn_els_on_nds(mn_mod.els, bdry_nds2);
els_in = tmp_els1 & ~tmp_els2; %restrict to elements INSIDE bdry_nds1 (i.e. not connected to bdry_nds2);
nds_in = setdiff(fn_nds_on_els(mn_mod.els, els_in), bdry_nds1);
all_els_in_subdomain_in_val = fn_all_els_inside_bdry(mn_mod.els, nds_in, els_in);
[~, ~, mn_mod.els, mn_mod.el_mat_i, mn_mod.el_abs_i, mn_mod.el_typ_i] = fn_remove_unused_elements(~all_els_in_subdomain_in_val, mn_mod.els, mn_mod.el_mat_i, mn_mod.el_abs_i, mn_mod.el_typ_i);
end