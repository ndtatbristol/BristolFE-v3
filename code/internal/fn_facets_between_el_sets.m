function [facet_nds, facet_cent, facet_area, facet_norm] = fn_facets_between_el_sets(nds, els,  el_typ_i, el_set_i, el_types, el_set_i_in, el_set_i_out)
%finds commone facets between el_set_i_in and el_set_i_out and returns
%details, with surface normals pointing out

%all faces on all elements in model (for unused nodes)
n_els = size(els, 1);
el_faces = fn_faces_from_els(els, (1:n_els)', el_typ_i, el_types);
max_faces = 0;
n_faces = 0;
for i = 1:numel(el_faces)
    max_faces = max(max_faces, size(el_faces{i}.fcs, 2));
    n_faces = n_faces + size(el_faces{i}.fcs, 1);
end
faces = zeros(n_faces, max_faces);
el_i = zeros(n_faces, 1);
i1 = 1;
for i = 1:numel(el_faces)
    i2 = i1 + size(el_faces{i}.fcs, 1) - 1;
    faces(i1:i2, 1:size(el_faces{i}.fcs, 2)) = el_faces{i}.fcs;
    el_i(i1:i2) = el_faces{i}.el_i;
    i1 = i1 + i2 + 1;
end

%reduce to just the faces in the two el_sets of interest
tmp = el_set_i(el_i);

i = tmp == el_set_i_in;
in_faces = faces(i, :);
in_el_i = el_i(i);

i = tmp == el_set_i_out;
out_faces = faces(i, :);
out_el_i = el_i(i);

in_faces_sorted = sort(in_faces, 2);
out_faces_sorted = sort(out_faces, 2);
[~, in_face_i, ~] = intersect(in_faces_sorted, out_faces_sorted, 'rows');

%finally the faces on the inside set on the boundary which should be
%consistently ordered because the order came from the element definitions
facet_nds = in_faces(in_face_i, :);

facet_cent = 0;
n = 0;
for i = 1:size(facet_nds, 2)
    n = n + double(facet_nds(:, i) > 0);
    facet_cent = facet_cent + nds(facet_nds(:, i), :);
end
facet_cent = facet_cent ./ n;
switch size(nds, 2)
    case 2
        tmp = nds(facet_nds(:, 1), :) - nds(facet_nds(:, 2), :); 
        facet_area = sqrt(sum(tmp .^ 2, 2));
        facet_norm = [-tmp(:, 2), tmp(:, 1)] ./ facet_area;
    case 3
        %two edge vectors
        d1 = nds(facet_nds(:, 2), :) - nds(facet_nds(:, 1), :);
        d2 = nds(facet_nds(:, 3), :) - nds(facet_nds(:, 1), :);

        facet_norm = cross(d1, d2, 2);
        facet_area = vecnorm(facet_norm, 2, 2);
        facet_norm = facet_norm ./ facet_area;
end

end