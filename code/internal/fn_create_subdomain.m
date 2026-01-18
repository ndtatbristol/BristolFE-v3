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

%Get elements in region - these guaranteed to be inside boundary and all
%boundary node layers
switch ndims
    case 2
        els_in_region = fn_2d_find_elements_in_region(dm_mod, inner_bdry_vtcs);
    case 3
        els_in_region = fn_3d_find_elements_in_region(dm_mod, inner_bdry_vtcs, inner_bdry_fcs);
end

method = 'new';
[dm_mod.bdry_lyrs, el_used] = fn_subdomain_bdry_layers(dm_mod.nds, dm_mod.els, els_in_region);

%Stick absorbing layer on
abs_layer_start_bdry_nds = dm_mod.nds(dm_mod.bdry_lyrs == 4, :);
el_centres = fn_calc_element_centres(dm_mod.nds, dm_mod.els);
%Restrict search to elements within possible range region + abs_layer_thick
cand_els = ~el_used & ...
    (all((el_centres < (max(abs_layer_start_bdry_nds) + abs_layer_thick)) & ...
    (el_centres > (min(abs_layer_start_bdry_nds) - abs_layer_thick)), 2));
dm_mod.el_abs_i(cand_els) = fn_quick_dist_to_point_bdry(el_centres(cand_els, :), abs_layer_start_bdry_nds) / abs_layer_thick;

els_in_use = el_used | cand_els;
els_in_use(dm_mod.el_abs_i > 1) = 0;

[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i] = fn_remove_unused_elements(els_in_use, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i);
[dm_mod.nds, dm_mod.els, old_nds, ~] = fn_remove_unused_nodes(dm_mod.nds, dm_mod.els);
dm_mod.main_nd_i = old_nds;
dm_mod.bdry_lyrs = dm_mod.bdry_lyrs(old_nds);


end

%--------------------------------------------------------------------------

function nds = fn_nds_on_els(els, els_to_consider)
nds = unique(els(els_to_consider, :));
end

function els = fn_els_on_nds(els, nds_to_consider)
els = any(ismember(els, nds_to_consider), 2);
end

function [nds_in, nds_on, nds_out, els_in, els_out] = fn_el_bdry_nds_and_els(els, els_in_bdry)
els_not_in_bdry =  ~els_in_bdry;
nds_on = intersect(fn_nds_on_els(els, els_in_bdry), fn_nds_on_els(els, els_not_in_bdry));
els_in = fn_els_on_nds(els, nds_on) & els_in_bdry;
els_out = fn_els_on_nds(els, nds_on) & els_not_in_bdry;
nds_in = setdiff(fn_nds_on_els(els, els_in), nds_on);
nds_out = setdiff(fn_nds_on_els(els, els_out), nds_on);
end

function bdry_fcs = fn_bdry_fcs(els, bdry_nds, el_typ_i, el_types)
%Inputs
%   els - n_els x n_max_nds_per_el matrix of elements to consider (it only
%   needs to be ones at bounday)
%   bdry_nds - vector of boundary nodes
%Returns
%   bdry_fcs - boundary face matrix indexed into bdry_nds?
el_i = (1:size(els, 1))';
el_faces = fn_faces_from_els(els, el_i, el_typ_i, el_types);
%first pass just to get max size of matrix to store results
max_nds_per_face = 0;
max_faces = 0;
for i = 1:numel(el_faces)
    max_nds_per_face = max(max_nds_per_face, size(el_faces{i}.fcs,2));
    max_faces = max_faces + size(el_faces{i}.fcs,1);
end
bdry_fcs = zeros(max_faces, max_nds_per_face);
%second pass - extract faces which only contain common nodes
k = 1;
for i = 1:numel(el_faces)
    j = all(ismember(el_faces{i}.fcs, bdry_nds), 2);
    bdry_fcs(k:k + nnz(j) - 1, 1:size(el_faces{i}.fcs, 2)) = el_faces{i}.fcs(j,:);
    k = k + nnz(j);
end
bdry_fcs(k:end, :) = [];
bdry_fcs = fn_unique_fcs(bdry_fcs);

[tf, idx] = ismember(bdry_fcs(:), bdry_nds);   % idx are positions in v for each element of m (linearized)
assert(all(tf), 'Some entries of m are not present in v.');
bdry_fcs = reshape(idx, size(bdry_fcs));   % same size as m, with indices into v
end

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
common_fcs = fn_unique_fcs(common_fcs);

[tf, idx] = ismember(common_fcs(:), common_nds);   % idx are positions in v for each element of m (linearized)
assert(all(tf), 'Some entries of m are not present in v.');

common_fcs = reshape(idx, size(common_fcs));   % same size as m, with indices into v


end
