function [bdry_lyrs, el_set_i] = fn_subdomain_bdry_layers_v2(nds, els, els_in_region, no_layers)

if~exist('no_layers', 'var')
    no_layers = 4;
end


bdry_lyrs = zeros(size(nds, 1), 1);
el_set_i = els_in_region;

%Move out one layer so that boundary (n1) will be on edge of els_in_region
[~, ~, ~, ~, e_out] = fn_el_bdry_nds_and_els(els, el_set_i);
el_set_i = el_set_i + e_out * 2;

for lyr = 1:no_layers
    %Work out and assign bdry nodes to layers
    [n1, n2_1, n3_1, ~, e2] = fn_el_bdry_nds_and_els(els, el_set_i ~= 0);
    bdry_lyrs(n1) = lyr; %first layer
    tmp = min(lyr + 2, no_layers + 1);
    el_set_i = el_set_i + e2 * tmp;
end

el_set_i(el_set_i == 0) = no_layers + 1;

% [n2_2, n3_2, n4, ~, e3] = fn_el_bdry_nds_and_els(els, el_used);
% n2 = union(n2_1, n2_2);
% bdry_lyrs(n2) = 2; %second layer
% n3 = union(n3_1, n3_2);
% bdry_lyrs(n3) = 3; %third layer
% bdry_lyrs(n4) = 4; %fourth layer
% el_used = el_used | e3;


%Order of elset_i for 4 bdry layers is
%[elset_i = 1] n1 [elset_i = 2] n2 [elset_i = 3] n3 [elset_i = 4] n4 [elset_i = 5]

end

