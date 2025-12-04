function dm_mod = fn_create_subdomain(mn_mod, el_types, inner_bdry_vtcs, inner_bdry_fcs, abs_layer_thick)
%Core function used for 2D and 3D. Only difference should be in
%determinining the interior elements.
%USAGE - 2D
%   dm_mod = fn_create_subdomain(mn_mod, inner_bdry_vtcs, [], abs_layer_thick)
%USAGE - 3D
%   dm_mod = fn_create_subdomain(mn_mod, inner_bdry_vtcs, inner_bdry_fcs, abs_layer_thick)

ndims = size(mn_mod.nds, 2);

%New version - based on tidier way of getting layers
%Needs to return something like this
    %             nds: [11424×2 double]
    %             els: [22516×3 double]
    %        el_mat_i: [22516×1 double]
    %        el_abs_i: [22516×1 double]
    %        el_typ_i: {22516×1 cell}
    %       bdry_lyrs: [11424×1 double]
    %       main_nd_i: [11424×1 double]
    % outer_bndry_pts: [882×2 double]
    % inner_bndry_pts: [361×2 double]

%First make a copy of the key parts of the main model
dm_mod = mn_mod;
if isfield(dm_mod, 'max_safe_time_step')
    dm_mod = rmfield(dm_mod, 'max_safe_time_step');
end
if isfield(dm_mod, 'design_centre_freq')
    dm_mod = rmfield(dm_mod, 'design_centre_freq');
end

%Create vector that will hold indices associating nodes with the 4 boundary
%layers in the subdomain
dm_mod.bdry_lyrs = zeros(size(mn_mod.nds, 1), 1);
% dm_mod.inner_bdry_vtcs = inner_bdry_vtcs;

%Get elements in region
switch ndims
    case 2
        el_used = fn_2d_find_elements_in_region(dm_mod, inner_bdry_vtcs);
    case 3
        el_used = fn_3d_find_elements_in_region(dm_mod, inner_bdry_vtcs, inner_bdry_fcs);
end

%Work out and assign bdry nodes to layers
for i = 1:4
    if i == 3 %for 3rd one, use boundary as start of absorbing region
        [el_i, common_nds, common_fcs] = fn_find_adjacent_els_to_els(dm_mod.els, mn_mod.el_typ_i, el_types, el_used, ~el_used);
        abs_layer_start_bdry_nds = dm_mod.nds(common_nds, :);
        abs_layer_start_bdry_fcs = common_fcs;
    else
        [el_i, common_nds] = fn_find_adjacent_els_to_els(dm_mod.els, mn_mod.el_typ_i, el_types, el_used, ~el_used);
    end
    dm_mod.bdry_lyrs(common_nds) = i;
    el_used(el_i) = 1;
end

% figure;c = 'rgbm'; h = fn_show_geometry(dm_mod, [], el_types, []);for i = 1:4; fn_plot_line(dm_mod.nds(dm_mod.bdry_lyrs == i,:), [c(i), '.']); end


%Delete original interface elements and then regenerate later in this 
%function - this is to avoid potential instability at edge of domain where 
%original interface elements copied from main mesh may now be on free edges
%Update: this doesn't seem to be necessary anymore now that instabilities fixed in general
% dm_mod = fn_remove_fluid_solid_interface_els(dm_mod);


%Add the absorbing layers by working out from centre of region
cand_els = ~el_used;
switch ndims
    case 2
        % dm_mod.el_abs_i(cand_els) = fn_dist_point_to_bdry_2D(fn_calc_element_centres(dm_mod.nds, dm_mod.els(cand_els, :)), abs_layer_start_bdry) / abs_layer_thick;
        dm_mod.el_abs_i(cand_els) = fn_2d_signed_dist_to_bdry(fn_calc_element_centres(dm_mod.nds, dm_mod.els(cand_els, :)), abs_layer_start_bdry_nds, abs_layer_start_bdry_fcs) / abs_layer_thick;
    case 3
        dm_mod.el_abs_i(cand_els) = fn_dist_point_to_bdry_3D(fn_calc_element_centres(dm_mod.nds, dm_mod.els(cand_els, :)), abs_layer_start_bdry_nds, abs_layer_start_bdry_fcs) / abs_layer_thick;
end
els_in_use = ones(size(dm_mod.els, 1), 1);
els_in_use(dm_mod.el_abs_i > 1) = 0;

[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i] = fn_remove_unused_elements(els_in_use, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i);
[dm_mod.nds, dm_mod.els, old_nds, ~] = fn_remove_unused_nodes(dm_mod.nds, dm_mod.els);
dm_mod.main_nd_i = old_nds;
dm_mod.bdry_lyrs = dm_mod.bdry_lyrs(old_nds);

%Reinstate fluid-solid interface elements - no longer needed if they are
%not deleted at start of process
% dm_mod = fn_add_fluid_solid_interface_els(dm_mod);

% free_ed = fn_find_free_edges(dm_mod.els);
% 
% dm_mod.outer_bndry_pts = [dm_mod.nds(free_ed, 1), dm_mod.nds(free_ed, 2)];
% dm_mod.int_el_i = fn_elements_in_region(dm_mod, dm_mod.inner_bndry_pts);

end

%--------------------------------------------------------------------------

function [adj_els, common_nds, common_fcs] = fn_find_adjacent_els_to_els(els, el_typ_i, el_types, els_to_consider, els_to_choose_from)
%returns logical array [size(els,1)x1] of elements in
%els(els_to_choose_from,:) that share common_nds with
%els(els_to_consider,:)

%unique nodes in els_to_consider
un_nds_to_consider = unique(els(els_to_consider,:));
un_nds_to_consider(un_nds_to_consider == 0) = [];

%unique nodes in els_to_choose_from
un_nds_to_choose_from = unique(els(els_to_choose_from,:));
un_nds_to_choose_from(un_nds_to_choose_from == 0) = [];

%nodes defining the boundary are those in both sets
common_nds = intersect(un_nds_to_consider, un_nds_to_choose_from);

%find which of the els_to_choose_from these are in
els_with_common_nodes = ismember(els, common_nds);
els_with_common_nodes(~els_to_choose_from, :) = 0;
adj_els = sum(els_with_common_nodes, 2) > 0;

if nargout < 3
    %return here if boundary faces are not required to save time
    return
end

%find list of unique faces that ONLY involve bdry nodes as these will
%describe boundary (in 2D and 3D)
el_i = (1:size(els, 1))';
el_faces = fn_faces_from_els(els, el_i, el_typ_i, el_types);
%first pass just to get max size of matrix to store results
max_nds_per_face = 0;
max_faces = 0;
for i = 1:numel(el_faces)
    max_nds_per_face = max(max_nds_per_face, size(el_faces{i}.fcs,2));
    max_faces = max_faces + size(el_faces{i}.fcs,1);
end
common_fcs = zeros(max_faces, max_nds_per_face);
%second pass - extract faces which only contain common nodes
k = 1;
for i = 1:numel(el_faces)
    j = all(ismember(el_faces{i}.fcs, common_nds), 2);
    common_fcs(k:k + nnz(j) - 1, 1:size(el_faces{i}.fcs, 2)) = el_faces{i}.fcs(j,:);
    k = k + nnz(j);
end
common_fcs(k:end, :) = [];
common_fcs = sort(common_fcs, 2);
common_fcs = unique(common_fcs, 'rows');


[tf, idx] = ismember(common_fcs(:), common_nds);   % idx are positions in v for each element of m (linearized)
assert(all(tf), 'Some entries of m are not present in v.');

common_fcs = reshape(idx, size(common_fcs));   % same size as m, with indices into v


end
