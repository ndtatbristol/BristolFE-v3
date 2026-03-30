function el_edges = fn_edges_from_els(els, el_i, el_typ_i, el_types)
%this returns a cell array (one cell per element type), each cell
%containing a matrix of element edge nodes for elements of that type in the
%model. Structure will look like:
%   el_faces{q}.el_typ = name of element given by el_types{q}
%   el_faces{q}.eds = matrix of nodes of edges of each element of el_types{q} with 
%       size [n_els_of_typ_q * n_eds_per_el_of_type_q, n_dim - 1]
%   el_faces{q}.el_i = vector of element numbers associated with each row
%       of el_faces{q}.eds of size [n_els_of_typ_q * n_eds_per_el_of_type_q, 1]
%For 2D elements this function returns faces (identical to fn_faces_from_els)

el_type_info = fn_el_type_info();
un_el_typ_i = unique(el_typ_i);
el_i = el_i(:);

for i = 1:numel(un_el_typ_i)
    ei = el_types{un_el_typ_i(i)};

    switch fn_el_shape(el_type_info, ei)
        case 'triangular' %2D triangular
            fc_i = [
                1, 2
                2, 3
                3, 1];
        case 'quadrilateral' %2D quadrilateral
            fc_i = [
                1, 2
                2, 3
                3, 4
                4, 1];
        case 'line' %2D line - note that for lines this function returns faces (same as fn_faceas_from_els)
            fc_i = [
                1, 2];
        case 'tetrahedral' %3D tetrahedral
            fc_i = [
                1, 2
                1, 3
                1, 4
                2, 3
                2, 4
                3, 4];
        case 'hexahedral' %3D hexahedral
            fc_i = [
                1, 2
                2, 3
                3, 4
                4, 1
                1, 5
                2, 6
                3, 7
                4, 8
                5, 6
                6, 7
                7, 8
                8, 5
                ];

    end
    j = el_typ_i == un_el_typ_i(i);

    el_edges{i}.el_typ_i = ei;
    el_edges{i}.eds = reshape(els(j, fc_i')', size(fc_i, 2), [])';
    el_edges{i}.el_i = reshape(reshape(el_i(j), 1, []) .* ones(size(fc_i, 1), 1), [], 1);
end

end
