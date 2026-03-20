function el_typ_i = fn_el_types_of_state(el_types, state)
%Returns indices of el_types that match specified state (e.g. solid, fluid
%etc)
states = cell(1,numel(el_types));
for i = 1:numel(el_types)
    tmp = fn_query_el_type_info(el_types{i});
    states{i} = tmp.state;
end
el_typ_i = find(strcmp(states, state));
end