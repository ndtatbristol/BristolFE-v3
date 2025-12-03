function shape = fn_el_shape(el_type_info, name)
%Returns the shape of the named element (triangular, hexahedral etc)
for i = 1:numel(el_type_info)
    if strcmp(el_type_info{i}.name, name)
        shape = el_type_info{i}.shape;
        return
    end
end
error(['Element ', name, ' not found'])

end