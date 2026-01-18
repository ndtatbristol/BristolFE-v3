function [bdry_lyrs, el_used] = fn_subdomain_bdry_layers(nds, els, els_in_region)

bdry_lyrs = zeros(size(nds, 1), 1);
el_used = els_in_region;



%Move out one layer so that boundary (n1) will be on edge of els_in_region
[~, ~, ~, ~, e_out] = fn_el_bdry_nds_and_els(els, el_used);
el_used = el_used | e_out;
%Work out and assign bdry nodes to layers
[n1, n2_1, n3_1, ~, e2] = fn_el_bdry_nds_and_els(els, el_used);
bdry_lyrs(n1) = 1; %first layer
el_used = el_used | e2;
[n2_2, n3_2, n4, ~, e3] = fn_el_bdry_nds_and_els(els, el_used);
n2 = union(n2_1, n2_2);
bdry_lyrs(n2) = 2; %second layer
n3 = union(n3_1, n3_2);
bdry_lyrs(n3) = 3; %third layer
bdry_lyrs(n4) = 4; %fourth layer
el_used = el_used | e3;


%rRder is
%[els_in_region] n1 [e1] n2 [e2] n3 [e3] n4 [absorbing elements]
%[el_used .............................] n4 [absorbing elements]

end

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

% function bdry_fcs = fn_bdry_fcs(els, bdry_nds, el_typ_i, el_types)
% %Inputs
% %   els - n_els x n_max_nds_per_el matrix of elements to consider (it only
% %   needs to be ones at bounday)
% %   bdry_nds - vector of boundary nodes
% %Returns
% %   bdry_fcs - boundary face matrix indexed into bdry_nds?
% el_i = (1:size(els, 1))';
% el_faces = fn_faces_from_els(els, el_i, el_typ_i, el_types);
% %first pass just to get max size of matrix to store results
% max_nds_per_face = 0;
% max_faces = 0;
% for i = 1:numel(el_faces)
%     max_nds_per_face = max(max_nds_per_face, size(el_faces{i}.fcs,2));
%     max_faces = max_faces + size(el_faces{i}.fcs,1);
% end
% bdry_fcs = zeros(max_faces, max_nds_per_face);
% %second pass - extract faces which only contain common nodes
% k = 1;
% for i = 1:numel(el_faces)
%     j = all(ismember(el_faces{i}.fcs, bdry_nds), 2);
%     bdry_fcs(k:k + nnz(j) - 1, 1:size(el_faces{i}.fcs, 2)) = el_faces{i}.fcs(j,:);
%     k = k + nnz(j);
% end
% bdry_fcs(k:end, :) = [];
% bdry_fcs = fn_unique_fcs(bdry_fcs);
%
% [tf, idx] = ismember(bdry_fcs(:), bdry_nds);   % idx are positions in v for each element of m (linearized)
% assert(all(tf), 'Some entries of m are not present in v.');
% bdry_fcs = reshape(idx, size(bdry_fcs));   % same size as m, with indices into v
% end
%
% %For OLD version
%
% function [adj_els, common_nds, common_fcs] = fn_find_adjacent_els_to_els(els, el_typ_i, el_types, els_to_consider, els_to_choose_from)
% %returns logical array [size(els,1)x1] of elements in
% %els(els_to_choose_from,:) that share common_nds with
% %els(els_to_consider,:)
%
% %unique nodes in els_to_consider
% un_nds_to_consider = unique(els(els_to_consider,:));
% un_nds_to_consider(un_nds_to_consider == 0) = [];
%
% %unique nodes in els_to_choose_from
% un_nds_to_choose_from = unique(els(els_to_choose_from,:));
% un_nds_to_choose_from(un_nds_to_choose_from == 0) = [];
%
% %nodes defining the boundary are those in both sets
% common_nds = intersect(un_nds_to_consider, un_nds_to_choose_from);
%
% %find which of the els_to_choose_from these are in
% els_with_common_nodes = ismember(els, common_nds);
% els_with_common_nodes(~els_to_choose_from, :) = 0;
% adj_els = sum(els_with_common_nodes, 2) > 0;
%
% if nargout < 3
%     %return here if boundary faces are not required to save time
%     return
% end
%
% %find list of unique faces that ONLY involve bdry nodes as these will
% %describe boundary (in 2D and 3D)
% el_i = (1:size(els, 1))';
% el_faces = fn_faces_from_els(els, el_i, el_typ_i, el_types);
% %first pass just to get max size of matrix to store results
% max_nds_per_face = 0;
% max_faces = 0;
% for i = 1:numel(el_faces)
%     max_nds_per_face = max(max_nds_per_face, size(el_faces{i}.fcs,2));
%     max_faces = max_faces + size(el_faces{i}.fcs,1);
% end
% common_fcs = zeros(max_faces, max_nds_per_face);
% %second pass - extract faces which only contain common nodes
% k = 1;
% for i = 1:numel(el_faces)
%     j = all(ismember(el_faces{i}.fcs, common_nds), 2);
%     common_fcs(k:k + nnz(j) - 1, 1:size(el_faces{i}.fcs, 2)) = el_faces{i}.fcs(j,:);
%     k = k + nnz(j);
% end
% common_fcs(k:end, :) = [];
% common_fcs = fn_unique_fcs(common_fcs);
%
% [tf, idx] = ismember(common_fcs(:), common_nds);   % idx are positions in v for each element of m (linearized)
% assert(all(tf), 'Some entries of m are not present in v.');
%
% common_fcs = reshape(idx, size(common_fcs));   % same size as m, with indices into v
%
%
% end