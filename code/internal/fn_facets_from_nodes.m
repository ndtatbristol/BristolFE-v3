function [facet_nds, facet_cent, facet_area, facet_norm] = fn_facets_from_nodes(nds, els, el_typ_i, el_types, nds_defining_facets)
%SUMMARY
%   Given connected set of nodes in mesh, this returns details of the
%   facets that they define. Each facet is defined by directly conncted nds
%   in nds_defining facets. Typically these may have come from a boundary
%   layer calculation for a subdomain model.

%potentially need to deal with (i.e. remove) interface elements here?

n_els = size(els,1);

el_faces = fn_faces_from_els(els, (1:n_els)', el_typ_i, el_types);

% local node numbering of tet faces

facet_nds = el_faces{1}.fcs;

q = fn_el_type_info(); face_def = q{strcmp(el_types, el_faces{1}.el_typ_i)}.faces;


% corresponding element numbers
face_els = repelem((1:n_els).',size(face_def,1));

% retain faces entirely on supplied surface
is_surf_face = all(ismember(facet_nds,nds_defining_facets),2);

facet_nds = facet_nds(is_surf_face,:);
face_els = face_els(is_surf_face);

% canonical ordering for comparison
facet_nds = sort(facet_nds, 2);

% find duplicates (shared by 2 elements)
[ufaces,~,ic] = unique(facet_nds, 'rows');
n_occ = accumarray(ic,1);

% keep only external faces
keep = n_occ(ic)==1;

facet_nds = facet_nds(keep,:);
face_els = face_els(keep);

%%%%%

% find duplicates (shared by 2 elements)
facet_nds = unique(facet_nds, 'rows');

switch size(nds, 2)
    case 2
        facet_nds = fn_2d_consistent_facet_nodes(facet_nds);
        facet_cent = (nds(facet_nds(:, 1), :) + nds(facet_nds(:, 2), :)) / 2;
        tmp = nds(facet_nds(:, 1), :) - nds(facet_nds(:, 2), :); 
        facet_area = sqrt(sum(tmp .^ 2, 2));
        facet_norm = [-tmp(:, 2), tmp(:, 1)] ./ facet_area;
    case 3
        facet_nds = fn_polygons_to_triangles(facet_nds);
        facet_nds = fn_3d_consistent_facet_nodes(facet_nds);
end


end