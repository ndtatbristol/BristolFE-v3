function mod = fn_remove_fluid_solid_interface_els(mod, el_types)
%SUMMARY
%   Removes interface elements from model. Probably does not need to be
%   called by  user. Is always called from fn_add_fluid_solid_interface_els
%   to remove any existing interface elements before adding new ones (otherwise
%   existing ones will be doubled up, which probably doubles coupling).

int_el_typ_i = fn_el_types_of_state(el_types, 'fluid_solid_interface');

if ~any(ismember(mod.el_typ_i, int_el_typ_i))
    %No interface elements defined, so none to remove!
    return
end

int_els = zeros(size(mod.els, 1), 1);
for i = 1:numel(int_el_typ_i)
    int_els(mod.el_typ_i == int_el_typ_i(i)) = 1;
end

els_in_use = ~int_els;
[~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);
end
