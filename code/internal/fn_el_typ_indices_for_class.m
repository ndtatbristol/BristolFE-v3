function int_el_typ_i = fn_el_typ_indices_for_class(el_types, el_class)
%Returns indices of element type indices corresponding to class

el_classes = fn_el_classes();
int_el_typ_i = [];
for i = 1:numel(el_classes.(el_class))
    int_el_typ_i = [int_el_typ_i, find(strcmp(el_types, el_classes.(el_class){i}))];
end
end

