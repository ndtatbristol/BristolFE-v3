function int_el_typ_i = fn_interface_element_type_indices(el_types)
%Returns indices of element sypes corresponding to interface elements

el_typ = fn_el_types();
int_el_typ_i = [];
for i = 1:numel(el_typ.fluid_solid_interface)
    int_el_typ_i = [int_el_typ_i, find(strcmp(el_types, el_typ.fluid_solid_interface{i}))];
end
end

