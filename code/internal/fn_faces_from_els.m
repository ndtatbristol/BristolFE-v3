function el_faces = fn_faces_from_els(els, el_i, el_typ_i, el_types)
%this returns a cell array (one cell per element type), each cell
%containing a matrix of element face nodes for elements of that type in the
%model. Structure will look like:
%   el_faces{q}.el_typ = name of element given by el_types{q}
%   el_faces{q}.fcs = matrix of nodes of faces of each element of el_types{q} with 
%       size [n_els_of_typ_q * n_fcs_per_el_of_type_q, n_dim - 1]
%   el_faces{q}.el_i = vector of element numbers associated with each row
%       of el_faces{q}.fcs of size [n_els_of_typ_q * n_fcs_per_el_of_type_q, 1]

el_type_info = fn_el_type_info();

un_el_typ_i = unique(el_typ_i);
el_i = el_i(:);

for i = 1:numel(un_el_typ_i)
    ei = el_types{un_el_typ_i(i)};

    switch fn_el_shape(el_type_info, ei)
        case 'triangular' %2D triangular
            fc_i = [
                1,2
                2,3
                3,1];
        case 'quadrilateral' %2d quadrilaterals
            fc_i = [
                1,2
                2,3
                3,4
                4,1];
        case 'line' %2D line
            fc_i = [
                1,2];
        case 'tetrahedral' %3D tetrahedral
            fc_i = [
                1,3,2
                1,2,4
                2,3,4
                1,4,3];
        case 'hexahedral' %3D hexahedral
            fc_i = [
                1,4,3,2
                1,2,6,5
                2,3,7,6
                3,4,8,7
                1,5,8,4
                5,6,7,8
                ];

    end
    j = el_typ_i == un_el_typ_i(i);

    el_faces{i}.el_typ_i = ei;
    el_faces{i}.fcs = reshape(els(j, fc_i')', size(fc_i, 2), [])';
    el_faces{i}.el_i = reshape(reshape(el_i(j), 1, []) .* ones(size(fc_i, 1), 1), [], 1);

end

end
